import sqlite3

from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app, send_from_directory,
    jsonify
)
from collections import defaultdict as dd
import datetime
import werkzeug
import secrets
from werkzeug.utils import secure_filename, escape
import email_validator
from email_validator import EmailNotValidError
from pscs.auth import login_required, is_logged_in
from pscs.db import get_db
from pscs.pscs import add_user_to_project, delete_data, check_user_permission
from werkzeug.security import check_password_hash, generate_password_hash
from pscs.messaging.mail import send_email
import pscs
import os
import string
from itsdangerous.url_safe import URLSafeTimedSerializer
from pscs.authtools.validation.registration import decode_token

bp = Blueprint("projects", __name__, url_prefix="/project")


@bp.route('/<id_project>', methods=['GET', 'POST'])
@login_required
def project(id_project):
    if request.method == 'GET':
        return display_private_project(id_project)
    elif request.method == 'POST':
        # Check for rename or delete
        if 'newName' in request.json:  # is rename
            new_project_name = request.json['newName']
            # Assert that user has permission
            has_perm = check_user_permission(permission_name='project_management',
                                             permission_value=1,
                                             id_project=id_project)
            if not has_perm:
                flash("You do not have permission to rename the project; contact the project manager.")
            elif has_perm:
                # Update name
                db = get_db()
                db.execute('UPDATE projects SET name_project = ? WHERE id_project = ?', (new_project_name, id_project))
                db.commit()
            return url_for('projects.project', id_project=id_project)
        elif 'delete' in request.json:  # is delete
            # Check that user is allowed
            has_perm = check_user_permission(permission_name='project_management',
                                             permission_value=1,
                                             id_project=id_project)
            if not has_perm:
                flash("You do not have permission to delete this project.")
                return url_for('projects.project', id_project=id_project)
            elif has_perm:
                db = get_db()
                pscs.db.delete_project(db, id_project=id_project)
        elif 'addUser' in request.json:
            # check that current user is allowed to add users
            has_perm = check_user_permission(permission_name='project_management',
                                             permission_value=1,
                                             id_project=id_project)
            if has_perm:
                # add user
                user_to_add = request.json['addUser']
                # get user id
                db = get_db()
                id_user = db.execute("SELECT id_user FROM users_auth WHERE name_user = ?", (user_to_add,)).fetchone()
                if id_user is None:
                    # user doesn't exist
                    flash(f"User {user_to_add} not found.")
                    return url_for('projects.project', id_project=id_project)
                add_user_to_project(id_user=id_user["id_user"], id_project=id_project, db=db,
                                    role='member', permissions={'data_read': 1})
            else:
                flash("You do not have permission to add users to the project.")
            return url_for("projects.project", id_project=id_project)
        elif "deleteData" in request.json:
            id_data = request.json["deleteData"]
            delete_data(id_data)
            return url_for("projects.project", id_project=id_project)
    return redirect(url_for('pscs.index'))


@bp.route('/<id_project>/public', methods=['GET', 'POST'])
def public_project(id_project):
    """
    Displays the public version of the project.
    Parameters
    ----------
    id_project : str
        ID of the project to display

    Returns
    -------

    """
    if request.method == "GET":
        # First check session if user has previously entered password
        if "project_review" in session.keys() and id_project in session["project_review"]:
            return display_public_project(id_project)  # password previously entered
        # Check whether the project is under peer review or is public
        db = get_db()
        public_status = db.execute("SELECT is_published, is_peer_review "
                                   "FROM projects "
                                   "WHERE id_project = ?", (id_project,)).fetchone()
        if public_status is None or (public_status["is_peer_review"] == 0 and public_status["is_published"] == 0):
            flash("The specified project either does not exist or is private.")
            return redirect(url_for("pscs.index"))
        if public_status["is_peer_review"]:
            # Prompt user for password
            return render_template("auth/prompt.html", prompt_label="Peer review password")
        elif public_status["is_published"]:
            return display_public_project(id_project)
        else:
            return redirect(url_for("pscs.index"))

    elif request.method == "POST":
        submitted_password = request.form["password"]
        db = get_db()
        # Get passhash from db
        project_passhash = db.execute("SELECT peer_password "
                                      "FROM projects_peer_review "
                                      "WHERE id_project = ?", (id_project,)).fetchone()
        if project_passhash is None:
            flash("There was an error with the project's password. Please contact PSCS admin.")
            return redirect(url_for("pscs.index"))
        project_passhash = project_passhash["peer_password"]
        if check_password_hash(project_passhash, submitted_password):
            # password is valid; display page
            append_to_session("project_review", id_project)
            return display_public_project(id_project)

        else:
            flash("Password is incorrect.")
            return redirect(url_for("pscs.public_project", id_project=id_project))
    return


def display_public_project(id_project):
    db = get_db()
    project_desc = db.execute("SELECT id_project, name_project, description "
                              "FROM projects "
                              "WHERE id_project = ? AND (is_peer_review=1 OR is_published=1)", (id_project,)).fetchone()
    project_summary = dict(project_desc)
    # get authors & affiliations
    project_summary["authors"] = db.execute("SELECT users_auth.name, users_auth.id_user, publication_authors.author_position "
                                            "FROM users_auth INNER JOIN publication_authors "
                                            "ON users_auth.id_user = publication_authors.id_user "
                                            "WHERE publication_authors.id_project = ? ", (id_project,)).fetchall()
    # Affiliations
    project_summary["affiliations"] = dict()
    for au in project_summary["authors"]:
        id_auth = au["id_user"]
        pos = au["author_position"]
        affil = db.execute("SELECT affiliation "
                           "FROM users_affiliation "
                           "WHERE id_user=?", (id_auth,)).fetchall()
        project_summary["affiliations"][pos] = ', '.join([aff['affiliation'] for aff in affil])

    # get external authors & affiliations
    external_authors = db.execute("SELECT name, email "
                                  "FROM external_author_info "
                                  "WHERE id_project = ?", (id_project,)).fetchall()
    # External authors are identified by their email + project
    for au in external_authors:
        ext_email = au["email"]
        # Get author info similar to internal author (name + author position)
        ext_pos = db.execute("SELECT author_position "
                             "FROM publication_external_authors "
                             "WHERE id_project = ? AND email = ?", (id_project, ext_email)).fetchone()["author_position"]
        ext_author_info = {"name": au["name"], "author_position": ext_pos}
        project_summary["authors"].append(ext_author_info)
        # Get affiliations
        external_affiliations = db.execute("SELECT affiliation "
                                           "FROM external_author_affiliation "
                                           "WHERE id_project = ? AND email = ?", (id_project, ext_email)).fetchall()
        project_summary["affiliations"][ext_pos] = ', '.join([aff["affiliation"] for aff in external_affiliations])

    # Sort author list by position
    project_summary["authors"].sort(key=lambda x: x["author_position"])
    # get analyses
    analyses = db.execute("SELECT id_analysis, analysis_name "
                          "FROM analysis "
                          "WHERE id_project = ? AND (is_peer_review=1 OR is_published=1)", (id_project,)).fetchall()
    project_summary["analyses"] = [dict(an) for an in analyses]
    public_analysis_id = set()
    for an in project_summary["analyses"]:
        public_analysis_id.add(an["id_analysis"])
    results_rows = db.execute("SELECT file_path, title, description, id_analysis "
                              "FROM results "
                              "WHERE id_project = ? AND (is_peer_review=1 OR is_published=1)", (id_project,)).fetchall()
    # Single DB fetch op means we have to sort them by analysis now
    project_summary["results"] = dd(list)
    for r in results_rows:
        an_id = r["id_analysis"]
        project_summary["results"][an_id].append(dict(r))
        project_summary["results"][an_id][-1]["file_name"] = os.path.basename(project_summary["results"][an_id][-1]["file_path"])

    return render_template("pscs/project_public.html",
                           project_summary=project_summary)


@bp.route('/<id_project>/results/<id_analysis>/<path:filename>', methods=['GET'])
def results(filename, id_project, id_analysis):
    # If logged in and has permission, don't need to check public
    if is_logged_in() and check_user_permission("data_read", 1, id_project):
        return private_results(filename, id_project, id_analysis)
    db = get_db()
    public_status = db.execute("SELECT is_published, is_peer_review "
                               "FROM projects "
                               "WHERE id_project = ?", (id_project,)).fetchone()
    if public_status is None:
        return  # Problem
    if public_status["is_published"]:
        return public_results(filename, id_project, id_analysis)
    elif public_status["is_peer_review"]:
        return review_results(filename, id_project, id_analysis)


# @login_required
def private_results(filename, id_project, id_analysis):
    if is_logged_in():
        has_perm = check_user_permission("data_read", 1, id_project)
        if not has_perm:
            return
        res_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis)
        return send_from_directory(res_dir, secure_filename(filename))
    else:
        return


def review_results(filename, id_project, id_analysis):
    if "project_review" in session.keys() and id_project in session["project_review"]:
        res_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis)
        return send_from_directory(res_dir, secure_filename(filename))
    return


def public_results(filename, id_project, id_analysis):
    res_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis)
    return send_from_directory(res_dir, secure_filename(filename))


def append_to_session(session_key: str,
                      value):
    """
    Appends `value` to the current session under the `session_key` list. If `session_key` is not defined, creates the
    list.
    Parameters
    ----------
    session_key : str
        Key for the list to which to append `value`.
    value
        Value to store.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If there is already a value under session_key and it is not a list.
    """
    if session_key not in session.keys():
        session[session_key] = []
    if not isinstance(session[session_key], list):
        raise ValueError(f"{session_key} is already used to store non-list data.")
    session[session_key].append(value)
    return


def display_private_project(id_project):
    """
    Renders the view for projects that have not yet been published.
    Parameters
    ----------
    id_project : str
        ID of the project to display

    Returns
    -------

    """
    db = get_db()
    user_read = check_user_permission('data_read', 1, id_project)
    user_write = check_user_permission('data_write', 1, id_project)
    id_user = g.user['id_user']
    if user_read:  # Check that user has read permission
        # Display files for this project only
        session['CURRENT_PROJECT'] = id_project
        project_name = db.execute('SELECT name_project FROM projects WHERE id_project = ?', (id_project,)).fetchone()
        project_name = project_name["name_project"]
        # Get analyses
        project_analyses = db.execute('SELECT id_analysis, analysis_name FROM analysis WHERE id_project = ?',
                                      (id_project,)).fetchall()
        analyses = {}  # used to store analysis name, keyed by ID
        analysis_nodes = {}
        for an in project_analyses:
            analyses[an['id_analysis']] = an['analysis_name']
            # Get input nodes with analysis
            project_inputs_db = db.execute('SELECT node_id, node_name FROM analysis_inputs WHERE id_analysis = ?',
                                           (an['id_analysis'],)).fetchall()
            project_inps = {}
            for inp in project_inputs_db:
                project_inps[inp['node_id']] = inp['node_name']
            analysis_nodes[an['id_analysis']] = project_inps
        # Get files associated with project
        project_data = db.execute('SELECT id_data, file_path FROM data WHERE id_project = ? AND id_user = ?',
                                  (id_project, id_user)).fetchall()
        files = {}
        for project_file in project_data:
            files[project_file['id_data']] = os.path.basename(project_file['file_path'])
        project_data_summary = db.execute(
            'SELECT id_data, file_path, data_type, file_hash, data_uploaded_time FROM data WHERE id_project = ?',
            (id_project,)).fetchall()

        # Get results
        results_rows = db.execute("SELECT file_path, title, description, id_analysis FROM results WHERE id_project = ?",
                                  (id_project,)).fetchall()
        results_files = []
        for r in results_rows:
            rr = dict(r)
            rr["file_name"] = os.path.basename(r["file_path"])
            results_files.append(rr)

        # Get users associated with project
        users = db.execute("SELECT name_user "
                           "FROM users_auth INNER JOIN projects_roles "
                           "ON users_auth.id_user = projects_roles.id_user WHERE id_project = ?",
                           (id_project,)).fetchall()
        summary = {"id_project": id_project, "id_user": id_user}
        user_list = []
        for u in users:
            user_list.append(u['name_user'])
        public_status = _get_project_publication_status(id_project)
        return render_template("pscs/project.html",
                               project_name=project_name,
                               analyses=analyses,
                               files=files,
                               analysis_nodes=analysis_nodes,
                               project_data_summary=project_data_summary,
                               results=results_files,
                               user_list=user_list,
                               project_status=public_status,
                               summary=summary)


@bp.route('/<id_project>/publish', methods=["GET", "POST"])
@login_required
def project_publish(id_project):
    if request.method == "GET":
        has_perm = check_user_permission(permission_name='project_management',
                                         permission_value=1,
                                         id_project=id_project)
        if not has_perm:
            flash("You do not have permission to do this; contact the project owner.")
            return redirect(url_for("projects.project", id_project=id_project))
        # Has permission to publish; display publication page
        db = get_db()
        project_info = db.execute("SELECT name_project, description "
                                  "FROM projects "
                                  "WHERE id_project = ?", (id_project,)).fetchone()
        # Author info
        authors, missing_name = _get_authors(id_project)

        # Data
        # List data associated with project
        project_data = _get_data(id_project)

        # Analyses & results
        analyses, result, unrun_analyses = _get_analyses(id_project)

        # Public status
        public_status = _get_project_publication_status(id_project)
        return render_template("pscs/project_publish.html",
                               proj=project_info,
                               authors=authors,
                               missing_name=missing_name,
                               project_data=project_data,
                               analyses=analyses,
                               results=result,
                               unrun_analyses=unrun_analyses,
                               public_status=public_status)
    elif request.method == "POST":
        # Before publishing, check that user is allowed, validate input, etc.
        has_perm = check_user_permission(permission_name="project_management",
                                         permission_value=1,
                                         id_project=id_project)
        if not has_perm:
            # Users will only end up here if their permissions have been changed while trying to publish a project
            flash("You do not have permission to publish this project.")
            return redirect(url_for("pscs.index"))

        # Verify that user has input the correct project name
        db = get_db()
        name_project = db.execute("SELECT name_project "
                                  "FROM projects "
                                  "WHERE id_project = ?", (id_project,)).fetchone()["name_project"]
        user_submmited_name = request.json["confirmation"]
        if name_project != user_submmited_name:
            flash("The project name you entered did not match the project's actual name. Please try again.")
            return redirect(url_for("projects.project_publish", id_project=id_project))

        # Checks have been passed; user is allowed and has entered the correct confirmation.
        # Validate the user input.
        publication_info = request.json
        # Now need to check all data
        # Four things are sent: publication_type, authorlist, analyses, data
        # Validate publication type
        valid_pubtype = validate_publication_type(publication_info["publication_type"], id_project=id_project)
        valid_authors, invalid_emails = validate_authorlist(publication_info["authorlist"], id_project=id_project)
        valid_analyses = validate_analyses(publication_info["analyses"], id_project=id_project)
        valid_data = validate_data(publication_info["data"], id_project=id_project)
        all_valid = valid_pubtype and valid_authors and valid_analyses and valid_data
        # Get author names
        if all_valid:
            # Fetch results associated with analyses
            publication_info["results"] = {}
            for an_id in publication_info["analyses"]:
                publication_info["results"][an_id] = db.execute("SELECT id_result "
                                                                "FROM results "
                                                                "WHERE id_analysis = ? AND id_project = ?", (an_id, id_project))
            if publication_info["publication_type"] == "peer review":
                # Log everything into db
                prepare_publication(id_project, publication_info)
                # Wait for external authors to approve before proceeding
                if _contains_external_authors(publication_info["authorlist"]):
                    _set_hold(id_project, commit=True)
                    # Notify external authors
                    notify_external_authors(id_project)
                    flash("The project is staged for peer review; it will be available for review once the external "
                          "authors confirm their association with the project.")
                else:
                    # Proceed with peer review
                    project_peer_review(id_project)
            elif publication_info["publication_type"] == "public":
                # Log everything into db
                prepare_publication(id_project, publication_info)
                # If project is currently in peer review, move forward. Otherwise, get confirmation from external author
                has_ext_auth = _contains_external_authors(publication_info["authorlist"])
                current_pub_status = _get_project_publication_status(id_project)
                if not has_ext_auth or (has_ext_auth and current_pub_status == "peer_review"):
                    project_public(id_project)
                else:
                    _set_hold(id_project, 1)
                    notify_external_authors(id_project)
                    flash("The project is staged to be public; it will be publicly viewable once the external "
                          "authors confirm their association with the project.")

            return_url = url_for('pscs.index')
            response_json = {"url": return_url, "publish_success": 1, "success_message": ""}
            return jsonify(response_json)
        else:
            # These only happen if the project changes while the user is trying to publish it, or sneaky user.
            flash_msg = ""
            if not valid_pubtype:
                flash_msg += f"Project can't be published as {publication_info['publication_type']} - "
            if not valid_authors:
                flash_msg += "Some of the authors are PSCS users that are not associated with the project. - "
            if not valid_analyses:
                flash_msg += "Some of the selected analyses can't be published. Ensure that they have been run first. - "
            if not valid_data:
                flash_msg += "Some of the selected data is not associated with the project. - "
            flash_msg = flash_msg[:-3]
            flash(flash_msg)
            return_url = url_for('pscs.index')
            response_json = {"url": return_url, "publish_success": 0, "success_message": flash_msg}
            return jsonify(response_json)
    return redirect(url_for("pscs.index"))


def _get_project_publication_status(id_project) -> str:
    """Returns the publication status of the specified project. Returns either "public", "peer review", or "private"."""
    db = get_db()
    pub_status = db.execute("SELECT is_published, is_peer_review, is_on_hold "
                            "FROM projects "
                            "WHERE id_project=?", (id_project,)).fetchone()
    if pub_status["is_published"] == 1:
        public_status = "public"
        if pub_status["is_on_hold"] == 1:
            public_status += " on_hold"
    elif pub_status["is_peer_review"] == 1:
        public_status = "peer review"
        if pub_status["is_on_hold"] == 1:
            public_status += " on_hold"
    else:
        public_status = "private"
    return public_status


def validate_publication_type(pubtype: str, id_project: str) -> bool:
    """Returns whether the project can be set to the specified publication type."""
    public_status = _get_project_publication_status(id_project)
    valid_pubtypes = ["private", "peer review", "public"]
    if pubtype == public_status:
        return False
    elif pubtype not in valid_pubtypes:
        return False
    elif public_status == "private":
        return True  # if it's private, any publication is fine
    elif public_status == "peer review" and pubtype == "public":
        return True  # if it's under peer review, only 'public' is valid
    else:
        # There is the rare circumstance where projects that are published or set for peer review need to be private
        # Users can do this only within the window specified in the data use agreement.
        if pubtype == "private":
            time_published = _get_project_publication_time(id_project)
            now_utc = int(datetime.datetime.now().astimezone(datetime.timezone.utc).strftime("%s"))
            return now_utc - time_published <= current_app.config["PUBLIC_RECANT_TIMEOUT"]
    return False  # default


def _get_project_publication_time(id_project: str) -> float:
    """Returns the epoch time for when a project was published."""
    db = get_db()
    time_published = db.execute("SELECT date_published "
                                "FROM projects "
                                "WHERE id_project = ?", (id_project,)).fetchone()["date_published"]
    return datetime.datetime.strptime(time_published, "%Y-%m-%d %H:%M:%S").timestamp()


def validate_authorlist(authorlist: list,
                        id_project: str) -> (bool, list):
    """
    Verifies that the input list is correct. Returns a list of incorrect external email addresses, if any.
    Parameters
    ----------
    authorlist : list
        List of dicts with either {"id_user": id} entries for PSCS users, or {"email": address} for external authors.
    id_project : str
        ID of the project to validate against.

    Returns
    -------
    bool
        Whether every entry in the list is valid.
    list
        List of invalid email addresses.
    """
    # Separate into PSCS and external
    pscs_ids = []
    external_emails = []
    for au in authorlist:
        try:
            email_validator.validate_email(au["id"])
            external_emails.append(au["id"])
        except EmailNotValidError:
            pscs_ids.append(au["id"])
    pscs_valid = _validate_pscs_authorlist(pscs_ids, id_project)
    external_valid, invalid_emails = _validate_external_authorlist(external_emails)
    return ((pscs_valid == 1) and (external_valid == 1), invalid_emails)


def _contains_external_authors(authorlist: list) -> bool:
    """Checks whether any of the authors in the list are external to PSCS."""
    for au in authorlist:
        if "email" in au.keys():
            return True
    return False


def _validate_external_authorlist(authorlist: list) -> (bool, list):
    """
    Checks whether the email addresses in the authorlist are valid.
    Parameters
    ----------
    authorlist : list
        List of email addresses belonging to external authors to contact.

    Returns
    -------
    bool
        Whether all entries in the list are valid.
    list
        List of invalid entries.
    """
    invalid_emails = []
    for au in authorlist:
        try:
            email_validator.validate_email(au)
        except EmailNotValidError:
            invalid_emails.append(au)
    return (len(invalid_emails) == 0, invalid_emails)


def _validate_pscs_authorlist(authorlist: list,
                              id_project: str) -> bool:
    """
    Checks whether the input PSCS users in the list match those that are associated with the project
    Parameters
    ----------
    authorlist : list
        List of user IDs for authors with PSCS accounts.
    id_project : str
        ID of the project that is being published.

    Returns
    -------
    bool
        Whether the list of authors is valid.
    """
    authorset = set(authorlist)
    # Get list of all authors associated with project
    db = get_db()
    authors = db.execute("SELECT id_user "
                         "FROM projects_roles "
                         "WHERE id_project = ?", (id_project,)).fetchall()
    author_ids = set([a["id_user"] for a in authors])
    # Input authorlist must contain all PSCS users. Input authorlist mustn't contain users unassociated with project
    return authorset == author_ids


def validate_analyses(analyses: list,
                      id_project: str) -> bool:
    """
    Verifies that the analyses are associated with the project.
    Parameters
    ----------
    analyses : list
        List of IDs of the analyses to be published.
    id_project : str
        ID of the project being published.

    Returns
    -------
    bool
        Whether the analyses are valid.
    """
    # Check that analyses are part of project
    analyses_set = set(analyses)
    db = get_db()
    analyses_db = db.execute("SELECT id_analysis "
                             "FROM analysis "
                             "WHERE id_project = ?", (id_project,)).fetchall()
    analyses_db_set = set([adb["id_analysis"] for adb in analyses_db])
    analyses_are_associated = analyses_set == analyses_set.intersection(analyses_db_set)

    # Check that analyses have been run
    results_db = db.execute("SELECT id_analysis "
                            "FROM results "
                            "WHERE id_project = ?", (id_project,)).fetchall()
    results_db_set = set([rdb["id_analysis"] for rdb in results_db])
    has_results = analyses_set.intersection(results_db_set) == analyses_set
    return analyses_are_associated and has_results


def validate_data(data: list,
                  id_project: str) -> bool:
    """
    Verifies that the data are associated with the project.
    Parameters
    ----------
    data : list
        List of IDs of the data to be published.
    id_project : str
        ID of the project being published.

    Returns
    -------
    bool
        Whether the list of data is valid.
    """
    data_set = set(data)
    db = get_db()
    data_db = db.execute("SELECT id_data "
                         "FROM data "
                         "WHERE id_project = ?", (id_project,)).fetchall()
    data_db_set = set([ddb["id_data"] for ddb in data_db])
    return data_set == data_set.intersection(data_db_set)


def _get_authors(id_project) -> (list, bool):
    """Returns the list of authors, each instance containing the name_user and name of the user. If a user doesn't have
    a name associated with their account, also returns True."""
    db = get_db()
    authors = db.execute("SELECT users_auth.name, users_auth.name_user, users_auth.id_user "
                         "FROM users_auth INNER JOIN projects_roles "
                         "ON users_auth.id_user = projects_roles.id_user "
                         "WHERE projects_roles.id_project = ?", (id_project,)).fetchall()
    author_dicts = {}
    for au in authors:
        author_dicts[au["id_user"]] = dict(au)

    # order the authors by author_position; this field is only available iff project was previously published
    pub_info = db.execute("SELECT id_user, author_position "
                          "FROM publication_authors "
                          "WHERE id_project = ?", (id_project,)).fetchall()
    if len(pub_info) > 0:
        for au_pub in pub_info:
            author_dicts[au_pub["id_user"]]["author_position"] = au_pub["author_position"]
    else:
        for id, _ in author_dicts.items():
            author_dicts[id]["author_position"] = 0
    author_list = [au for _, au in author_dicts.items()]

    # Get external authors, if any
    external_authors = db.execute("SELECT email, author_position "
                                  "FROM publication_external_authors "
                                  "WHERE id_project = ? AND confirmed = 1", (id_project,)).fetchall()
    if len(external_authors) > 0:
        for au in external_authors:
            au_name = db.execute("SELECT name "
                                 "FROM external_author_info "
                                 "WHERE id_project = ? AND email = ?", (id_project, au["email"])).fetchone()["name"]
            au_summ = {"name": au_name, "name_user": au["email"], "author_position": au["author_position"], "id_user": None}
            author_list.append(au_summ)

    author_list.sort(key=lambda x: x["author_position"])
    missing_name = False
    for a in authors:
        if a["name"] is None or len(a["name"]) == 0:
            missing_name = True
            break
    return author_list, missing_name


def _get_data(id_project) -> list:
    """Returns the list of data that has been uploaded to the project."""
    db = get_db()
    return db.execute("SELECT id_data, file_path, data_uploaded_time "
                      "FROM data "
                      "WHERE id_project = ?", (id_project,)).fetchall()


def _get_analyses(id_project):
    """Gets the analyses and results associated with each. If there are any analyses that were not run, return them
    separately."""
    db = get_db()
    # Get only analyses with results
    analyses_with_results = db.execute("SELECT analysis.id_analysis, analysis.analysis_name, results.id_result, results.file_path  "
                                       "FROM results INNER JOIN analysis "
                                       "ON results.id_analysis = analysis.id_analysis "
                                       "WHERE results.id_project = ?", (id_project,)).fetchall()
    analyses_summary = []
    results = dd(list)
    id_analyses_with_results = set()
    for awr in analyses_with_results:
        id_an = awr["id_analysis"]
        results[id_an].append({"id_result": awr["id_result"], "file_path": awr["file_path"]})
        if id_an in id_analyses_with_results:
            continue
        analyses_summary.append(awr)
        id_analyses_with_results.add(id_an)

    # Get analyses without results
    all_analyses = db.execute("SELECT id_analysis, analysis_name "
                              "FROM analysis "
                              "WHERE id_project = ?", (id_project,)).fetchall()
    analyses_without_results_summary = []
    for aa in all_analyses:
        id_an = aa["id_analysis"]
        if id_an in id_analyses_with_results:
            continue
        analyses_without_results_summary.append(aa)

    return analyses_summary, results, analyses_without_results_summary


def prepare_publication(id_project: str,
                        publication_info: dict):
    """
    Sets a project to peer review. Selected data, analyses, and results are placed behind a password prompt that does
    not require login.
    Parameters
    ----------
    id_project : str
        ID of the project to be published.
    publication_info : dict
        Dict containing what from the project should be made available.
    Returns
    -------
    None
    """
    db = get_db()
    if publication_info["publication_type"] == "peer review":
        peer, public = 1, 0
    elif publication_info["publication_type"] == "public":
        peer, public = 0, 1
    else:
        raise ValueError(f"Invalid publication type: {publication_info['publication_type']}")
    # Mark project as under peer review
    _set_project_status(id_project, is_peer_review=peer, is_published=public)
    # Mark each datum
    for id_data in publication_info["data"]:
        _set_data_status(id_data, is_peer_review=peer, is_published=public)
    # Mark each analysis
    for id_analysis in publication_info["analyses"]:
        _set_analysis_status(id_analysis, is_peer_review=peer, is_published=public)
    # Mark each result
    for id_an, results_list in publication_info["results"].items():
        for id_result in results_list:
            _set_result_status(id_result["id_result"], is_peer_review=peer, is_published=public)

    # Store authors in order
    pscs_authors, external_authors = _split_authlist(publication_info["authorlist"])
    _store_pscs_authors(id_project, pscs_authors, drop_previous=True)
    _store_external_authors(id_project, external_authors, drop_previous=True)

    # We've now done everything we need to do to prepare. Pass to next step.
    db.commit()
    return


def project_peer_review(id_project: str):
    """
    Publishes the specified project to be peer review. While under peer review, projects can be anonymously viewed with
    the supplied password.
    Parameters
    ----------
    id_project : str
        ID of the project to set to peer review.

    Returns
    -------
    None
    """
    # Remove hold
    db = get_db()
    _set_hold(id_project, value=0)

    peer_review_password = generate_password()
    # For placing password into db
    pass_hash = generate_password_hash(peer_review_password)
    # Check whether to update or insert; we're checking the table twice, but it's small
    previous = db.execute("SELECT id_project "
                          "FROM projects_peer_review "
                          "WHERE id_project=?", (id_project,)).fetchone()
    if previous is None:  # isn't already there; insert
        db.execute("INSERT INTO projects_peer_review "
                   "(id_project, peer_password) "
                   "VALUES(?,?)", (id_project, pass_hash))
    else:
        # Overwrite previous record
        db.execute("UPDATE projects_peer_review "
                   "SET peer_password = ? "
                   "WHERE id_project=?", (pass_hash, id_project))
    # Password is destined for HTML email; need to escape characters
    peer_review_password_html = werkzeug.utils.escape(peer_review_password)

    # Prepare email to notify user
    from pscs.templates.misc import peer_review_password_email
    project_info = db.execute("SELECT name_project, description "
                              "FROM projects "
                              "WHERE id_project = ?", (id_project,)).fetchone()
    project_url = "https://" + current_app.config["CURRENT_URL"] + url_for("projects.public_project", id_project=id_project)
    project_url = f"<a href='{project_url}'>{project_url}</a>"  # to make email content clickable
    peer_email = peer_review_password_email.format(title=project_info["name_project"],
                                                   description=project_info["description"],
                                                   project_url=project_url,
                                                   peer_review_password=peer_review_password_html)
    author_addresses = get_author_emails(id_project)
    send_email(author_addresses, subject="PSCS Project Set for Peer Review", body=peer_email)
    db.commit()
    return


def project_public(id_project):
    """
    Publishes the specified project to be publicly viewable.
    Parameters
    ----------
    id_project : str
        ID of the project to be made public.
    Returns
    -------
    None
    """
    db = get_db()
    # Remove hold, if any:
    _set_hold(id_project, value=0)

    # Remove peer review password, if any
    _remove_peer_review_password(id_project)
    from pscs.templates.misc import publication_notification
    project_info = get_project_summary(id_project)
    project_url = build_full_url(url_for("projects.public_project", id_project=id_project))
    project_url_html = f"<a href='{project_url}'>{project_url}</a>"
    notif_email = publication_notification.format(name_project=project_info["name_project"],
                                                  project_url=project_url_html)
    author_addresses = get_author_emails(id_project)
    send_email(author_addresses, subject="PSCS Project Publication", body=notif_email)
    db.commit()
    return


def _remove_peer_review_password(id_project, commit=False):
    """Removes the project's password associated with peer review, if any."""
    db = get_db()
    db.execute("DELETE FROM projects_peer_review WHERE id_project = ?", (id_project,))
    if commit:
        db.commit()
    return


def generate_password(length: int = 32) -> str:
    """Generates a password"""
    charset = string.digits + string.ascii_letters + string.punctuation
    return ''.join(secrets.choice(charset) for _ in range(length))


# NOTE: These are constructed as separate functions to avoid parameterizing the SQL string.
def _set_project_status(id_project, is_peer_review=None, is_published=None, commit=False):
    """Sets the publication status bits (is_peer_review, is_published) for the specified project. If None, their value
    is unchanged."""
    db = get_db()
    if is_peer_review is not None and is_published is not None:
        db.execute("UPDATE projects "
                   "SET is_peer_review = ?, is_published = ? "
                   "WHERE id_project = ?", (is_peer_review, is_published, id_project))
    elif is_peer_review is not None:
        db.execute("UPDATE projects SET is_peer_review = ? WHERE id_project = ?", (is_peer_review, id_project))
    elif is_published is not None:
        db.execute("UPDATE projects SET is_published = ? WHERE id_project = ?", (is_published, id_project))
    if commit:
        db.commit()
    return


def _set_data_status(id_data, is_peer_review=None, is_published=None, commit=False):
    """Sets the publication status bits (is_peer_review, is_published) for the specified datum. If None, their value is
    unchanged."""
    db = get_db()
    if is_peer_review is not None and is_published is not None:
        db.execute("UPDATE data "
                   "SET is_peer_review = ?, is_published = ? "
                   "WHERE id_data = ?", (is_peer_review, is_published, id_data))
    elif is_peer_review is not None:
        db.execute("UPDATE projects SET is_peer_review = ? WHERE id_data = ?", (is_peer_review, id_data))
    elif is_published is not None:
        db.execute("UPDATE projects SET is_published = ? WHERE id_data = ?", (is_published, id_data))
    if commit:
        db.commit()
    return


def _set_analysis_status(id_analysis, is_peer_review=None, is_published=None, commit=False):
    """Sets the publication bits (is_peer_review, is_published) for the specified analysis. If None, their value is
    unchanged."""
    db = get_db()
    if is_peer_review is not None and is_published is not None:
        db.execute("UPDATE analysis "
                   "SET is_peer_review = ?, is_published = ? "
                   "WHERE id_analysis = ?", (is_peer_review, is_published, id_analysis))
    elif is_peer_review is not None:
        db.execute("UPDATE analysis SET is_peer_review = ? WHERE id_analysis = ?", (is_peer_review, id_analysis))
    elif is_published is not None:
        db.execute("UPDATE analysis SET is_published = ? WHERE id_analysis = ?", (is_published, id_analysis))
    if commit:
        db.commit()
    return


def _set_result_status(id_result, is_peer_review=None, is_published=None, commit=False):
    """Sets the publication bits (is_peer_review, is_published) for the specified result. If None, the values are left
    unchanged."""
    db = get_db()
    if is_peer_review is not None and is_published is not None:
        db.execute("UPDATE results "
                   "SET is_peer_review = ?, is_published = ? "
                   "WHERE id_result = ?", (is_peer_review, is_published, id_result))
    elif is_peer_review is not None:
        db.execute("UPDATE results SET is_peer_review = ? WHERE id_result = ?", (is_peer_review, id_result))
    elif is_published is not None:
        db.execute("UPDATE results SET is_published = ? WHERE id_result = ?", (is_published, id_result))
    if commit:
        db.commit()
    return

def _set_result_peer_review(id_result, value=1, commit=False):
    """Sets is_peer_review to value in the DB for the analysis."""
    db = get_db()
    db.execute("UPDATE results SET is_peer_review = ? WHERE id_result = ?", (value, id_result))
    if commit:
        db.commit()


def _split_authlist(authlist: list) -> (list, list):
    """Splits the author list containing dicts into PSCS & external authors. The returned lists contain dicts with the
    "position" key indicating the order in which the author should appear."""
    pscs_authors = []
    ext_authors = []
    for position, au in enumerate(authlist):
        au["position"] = position
        id_is_email = False
        try:
            email_validator.validate_email(au["id"])
            id_is_email = True
        except EmailNotValidError:
            id_is_email = False
        if id_is_email:
            au["email"] = au["id"]
            ext_authors.append(au)
        else:
            pscs_authors.append(au)
    return pscs_authors, ext_authors


def _store_pscs_authors(id_project: str, pscs_authors: list, drop_previous: bool = True, commit=False):
    """
    Stroes the PSCS authors into the db for the specified project.
    Parameters
    ----------
    id_project : str
        ID of the project the authors are associated with.
    pscs_authors : list
        List of dicts with keys 'id' for the PSCS user id and 'position' for the author's position in the author list.
    drop_previous : bool
        Whether to drop previous authors associated with the project, if any.
    commit : bool
        Whether to commit changes to the DB.
    Returns
    -------
    None
    """
    db = get_db()
    if drop_previous:
        db.execute("DELETE FROM publication_authors WHERE id_project = ?", (id_project,))
    for au in pscs_authors:
        db.execute("INSERT INTO publication_authors "
                   "(id_project, id_user, author_position) "
                   "VALUES (?,?,?)", (id_project, au["id"], au["position"]))
    if commit:
        db.commit()
    return


def _store_external_authors(id_project: str, ext_authors: list, drop_previous: bool = True, commit: bool = False):
    """
    Stores info about external authors into the db for the specified project.
    Parameters
    ----------
    id_project : str
        ID of the project that the authors are associated with.
    ext_authors : list
        List of dicts with keys 'email' for the author's email and 'position' for the author's position in the author
        list.
    drop_previous : bool
        Whether to drop the previous authors associated with the project, if any.
    commit : bool
        Whether to commit changes to the DB.

    Returns
    -------
    None
    """
    db = get_db()
    if drop_previous:
        db.execute("DELETE FROM publication_external_authors WHERE id_project = ?", (id_project,))
    for au in ext_authors:
        db.execute("INSERT INTO publication_external_authors "
                   "(id_project, email, author_position) "
                   "VALUES (?,?,?)", (id_project, au["email"], au["position"]))
    if commit:
        db.commit()


def _set_hold(id_project: str, value=1, commit=False):
    """Notes that the project is being held (e.g., waiting for external users)."""
    db = get_db()
    db.execute("UPDATE projects SET is_on_hold = ? WHERE id_project = ?", (value, id_project))
    if commit:
        db.commit()


def get_author_emails(id_project: str) -> list:
    """Returns the emails of authors associated with the specified project."""
    db = get_db()
    # Get PSCS users
    pscs_authors = db.execute("SELECT users_auth.email "
                              "FROM users_auth INNER JOIN publication_authors "
                              "ON users_auth.id_user = publication_authors.id_user "
                              "WHERE publication_authors.id_project = ?", (id_project,)).fetchall()
    external_authors = db.execute("SELECT email "
                                  "FROM publication_external_authors "
                                  "WHERE id_project = ?", (id_project,)).fetchall()
    return [au["email"] for au in pscs_authors + external_authors]


def notify_external_authors(id_project: str):
    """Notifies external authors that they should input their information."""
    db = get_db()
    external_authors = db.execute("SELECT email "
                                  "FROM publication_external_authors "
                                  "WHERE id_project = ?", (id_project,)).fetchall()
    addresses = [au["email"] for au in external_authors]
    from pscs.templates.misc import external_author_info_request
    for addr in addresses:
        url_signer = URLSafeTimedSerializer(secret_key=current_app.config["SECRET_KEY"], salt="external_author")
        token = url_signer.dumps(addr)
        author_url = url_for("projects.external_author_info", id_project=id_project, token=token)
        external_author_url = build_full_url(author_url)
        notif_email = external_author_info_request.format(pscs_url=current_app.config["CURRENT_URL"],
                                                          author_info_url=external_author_url)
        subject = "PSCS Project Publication - Authorship Information Required"
        send_email([addr], subject=subject, body=notif_email)
    return


def build_full_url(local_url: str) -> str:
    """Builds the full URL including https://[domain] before the local URL."""
    return "https://" + current_app.config["CURRENT_URL"] + local_url


@bp.route("<id_project>/external_author/<token>", methods=["GET", "POST"])
def external_author_info(id_project, token):
    """
    Page for requesting an external author's information for publishing a project.
    Parameters
    ----------
    id_project : str
        ID of the project being published.
    token : str
        Signed token containing the external author's email.

    Returns
    -------

    """
    # Decode token to identify user
    url_signer = URLSafeTimedSerializer(secret_key=current_app.config["SECRET_KEY"], salt="external_author")
    external_author_email = url_signer.loads(token, max_age=current_app.config["EXTERNAL_AUTHOR_TIMEOUT_SECONDS"])
    # See if info matches pending project
    db = get_db()
    proj = db.execute("SELECT id_project "
                      "FROM publication_external_authors "
                      "WHERE id_project = ? AND email = ?", (id_project, external_author_email)).fetchone()
    if proj is None:  # id_project & email don't match user
        flash("The URL is invalid.")
        return redirect(url_for("pscs.index"))
    if request.method == "GET":
        # Get project info to display to user
        project_summary = get_project_summary(id_project)
        # Request user info
        return render_template("auth/external_author.html", project_summary=project_summary)
    elif request.method == "POST":
        # User has returned data
        data = request.form
        author_name = escape(data["name"])
        affiliations = data["affiliations"].splitlines()
        affiliations = [escape(aff) for aff in affiliations]
        noPHI = request.form["noPHI"] == "on"
        dataUse = request.form["dataUse"] == "on"
        flash_msg = ""
        # Form validation
        if not noPHI or not dataUse:
            flash_msg = "You must confirm that you have read, understood, and agreed to the PHI conditions and Data Use Agreement."
        if len(author_name) == 0:
            flash_msg = ' - '.join([flash_msg, "You must supply your name."])
        if len(affiliations) == 0:
            flash_msg = ' - '.join([flash_msg, "You must supply an affiliation. If unaffiliated, enter ""Unaffilliated."""])
        if len(flash_msg) != 0:
            flash(flash_msg)
            return redirect(url_for("projects.external_author_info", id_project=id_project, token=token))

        # Collect external author info
        author_info = {"email": external_author_email,
                       "name": author_name,
                       "affiliations": affiliations,
                       "noPHI": noPHI,
                       "dataUse": dataUse,
                       "ip": request.remote_addr}
        store_external_author_info(id_project, author_information=author_info)
        done = proceed_if_possible(id_project)
        if done:
            flash("Project is now ready for peer review.")
        else:
            flash("Project is now waiting for other external authors. An email notification will be sent when the "
                  "project is ready.")
        return redirect(url_for("pscs.index"))
    return


def proceed_if_possible(id_project: str):
    """Proceeds with publication stage if all external authors have approved."""
    # Check if all external authors have approved
    db = get_db()
    confirmed = db.execute("SELECT confirmed "
                           "FROM publication_external_authors "
                           "WHERE id_project = ?", (id_project,)).fetchall()
    # Checks if any of the authors have not confirmed.
    missing_confirmation = False in list(map(lambda x: x["confirmed"], confirmed))
    if missing_confirmation:
        return False  # not everyone has confirmed; exit
    # Move forward!
    pub_type = _get_project_publication_status(id_project)
    if "peer review" in pub_type:
        project_peer_review(id_project)
    elif "public" in pub_type:
        raise NotImplementedError("Setting projects to public not yet implemented.")
    return True


def store_external_author_info(id_project: str,
                               author_information: dict):
    """
    Stores the external author information into the DB.
    Parameters
    ----------
    id_project : str
        ID of the project for the author.
    author_information : dict
        Dict containing the following fields: email, name, affiliations, noPHI, dataUse, ip

    Returns
    -------
    None
    """
    db = get_db()
    db.execute("INSERT INTO external_author_info "
               "(id_project, email, name, confirmed_datause, confirmed_phi, ip) "
               "VALUES (?,?,?,?,?,?)", (id_project, author_information["email"],
                                        author_information["name"], author_information["dataUse"],
                                        author_information["noPHI"], author_information["ip"]))
    for aff in author_information["affiliations"]:
        db.execute("INSERT INTO external_author_affiliation "
                   "(id_project, email, affiliation) "
                   "VALUES (?,?,?)", (id_project, author_information["email"], aff))
    # Mark as confirmed
    db.execute("UPDATE publication_external_authors "
               "SET confirmed = 1 "
               "WHERE id_project = ? and email = ?", (id_project, author_information["email"]))
    db.commit()
    return


def get_project_summary(id_project: str) -> dict:
    """
    Builds a summary of the specified project.
    Parameters
    ----------
    id_project : str
        ID of the project to summarize.
    Returns
    -------
    dict
        Summary of the project.
    """
    db = get_db()
    proj_info = db.execute("SELECT name_project, description "
                           "FROM projects "
                           "WHERE id_project = ?", (id_project,)).fetchone()
    return proj_info

