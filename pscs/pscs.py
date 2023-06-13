import random

from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app, send_from_directory,
    jsonify
)
from markdown import markdown
from werkzeug.exceptions import abort
import werkzeug.utils
from werkzeug.utils import secure_filename
from werkzeug.security import check_password_hash, generate_password_hash
from pscs.auth import login_required, is_logged_in
import pscs.db
from pscs.db import get_db, get_unique_value_for_field
from pscs.metadata.metadata import get_metadata
from pscs.analysis.pipeline import node_parser
import os
import uuid
import numpy as np
from pscs.transfers.dispatching import dispatch, can_user_submit
from pscs.messaging.mail import send_email
from flask import send_from_directory
import plotly.express as px
import plotly
import pandas as pd
import json
import hashlib
import pathlib
import sqlite3
from collections import defaultdict as dd
import secrets
import string

bp = Blueprint("pscs", __name__)
ALLOWED_EXTENSIONS = {'csv', 'tsv'}
PATH_KEYWORD = 'path'


@bp.route('/')
def index():
    db = get_db()
    # Get project metadata to list for user
    posts = db.execute("SELECT title, body "
                               "FROM posts "
                               "ORDER BY date_created DESC "
                               "LIMIT 5").fetchall()
    posts_markdown = []
    for p in posts:
        post = {"title": p["title"]}
        post["body"] = markdown(p["body"])
        posts_markdown.append(post)
        # posts_markdown.append({p["title"]: markdown(p["body"])})
    return render_template("pscs/index.html", posts=posts_markdown)


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route(f"/upload/<id_project>/<name>", methods=['GET'])
@login_required
def download_file(id_project, name):
    return send_from_directory(current_app.config['UPLOAD_FOLDER'].format(id_project=id_project), name)


@bp.route('/upload', methods=['GET','POST'])
@login_required
def upload():
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            # Check that user has permission
            has_perm = check_user_permission(permission_name="data_write",
                                             permission_value=1,
                                             id_project=session["CURRENT_PROJECT"])
            if not has_perm:
                flash("You do not have permission to upload files.")
                return redirect(url_for("pscs.project", id_project=session["CURRENT_PROJECT"]))
            filename = secure_filename(file.filename)
            out_path = current_app.config['UPLOAD_FOLDER'].format(id_project=session["CURRENT_PROJECT"])
            os.makedirs(out_path, exist_ok=True)
            out_file = os.path.join(out_path, filename)
            file.save(out_file)
            if filename.endswith('tsv') or filename.endswith('csv'):
                meta_dict = get_metadata(out_file)
                sample_count, gene_count = meta_dict['table_dimensions']
                hash_value = meta_dict['table_hash']
                id_project = session['CURRENT_PROJECT']
                db = get_db()
                id_data = get_unique_value_for_field(db, "id_data", "data")
                db.execute(
                    'INSERT INTO data (id_data, id_user, id_project, file_path, data_type, file_hash)'
                    ' VALUES (?,?,?,?,?,?)',
                    (id_data, g.user['id_user'], id_project, out_file, 'table', hash_value))
                num_files = len(db.execute("SELECT id_data FROM data WHERE id_project = ?", (id_project,)).fetchall())
                db.execute(f'UPDATE projects SET num_files = {num_files} WHERE id_project = ?',
                           (id_project,))
                db.commit()
            return redirect(url_for('pscs.project', id_project=session['CURRENT_PROJECT']))
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <input type=file name=file>
      <input type=submit value=Upload>
    </form>
    '''


@bp.route('/projects_summary', methods=['GET'])
@login_required
def projects_summary():
    db = get_db()
    if g.user is not None:
        # Get project meta data to list for user
        user_projects = db.execute("SELECT projects.name_project, projects.description, projects.num_files, "
                                   "projects.id_project, projects.num_members, projects_roles.role "
                                   "FROM projects INNER JOIN projects_roles "
                                   "ON projects.id_project=projects_roles.id_project "
                                   "AND projects_roles.id_user=?", (g.user['id_user'],))
        return render_template('pscs/projects_summary.html', projects=user_projects)


@bp.route('/profile', methods=['GET'])
@login_required
def profile():
    if request.method == "GET":
        # Get basic user info
        id_user = g.user["id_user"]
        db = get_db()
        user_info = db.execute("SELECT name_user, email, creation_time_user, confirmed "
                               "FROM users_auth "
                               "WHERE id_user = ?", (id_user,)).fetchone()

        return render_template("pscs/profile.html", user_info=user_info)
    return render_template(url_for("pscs.index"))


@bp.route('/create_project', methods=['GET', 'POST'])
@login_required
def create_project():
    if request.method == 'POST':  # data has been sent to us
        project_name = request.form['name']
        project_description = request.form['description']
        error = None
        if project_name is None:
            error = "Project must have a name"
        if error is not None:
            flash(error)
        else:
            db = get_db()
            # Check whether user is allowed to make more
            can_create = below_project_limits(db, g.user["id_user"])
            if not can_create:
                flash("Project limits for user have been exceeded.")
                return redirect("create_project")
            # Insert into projects table
            id_project = get_unique_value_for_field(db, field="id_project", table="projects")

            # Create project
            db.execute("INSERT INTO projects (id_project, id_user, name_project, description) VALUES (?,?,?,?)",
                       (id_project, g.user['id_user'], project_name, project_description))
            # Add user to project
            add_user_to_project(id_user=g.user['id_user'], id_project=id_project, db=db, role='admin',
                                permissions={'data_read': 1, 'data_write': 1, 'project_management': 1,
                                             "analysis_read": 1, "analysis_write": 1, "analysis_execute": 1})
            db.commit()
            proj_dir = pathlib.Path(current_app.config['PROJECTS_DIRECTORY'].format(id_project=id_project))
            proj_dir.mkdir(exist_ok=True)
            results_dir = pathlib.Path(os.path.join(current_app.config['PROJECTS_DIRECTORY'], 'results').format(id_project=id_project))
            results_dir.mkdir(exist_ok=True)
            return redirect(url_for('pscs.project', id_project=id_project))
    return render_template("pscs/create.html")


def below_project_limits(db, id_user: str) -> bool:
    """
    Checks whether project limits have been exceeded by the user.
    Parameters
    ----------
    db : sqlite3.Connection
        Database connection to check.
    id_user : str
        Id of the user to check.

    Returns
    -------
    bool
        Whether the user is below project limits
    """
    project_count = len(db.execute("SELECT id_project "
                                   "FROM projects "
                                   "WHERE id_user = ?", (id_user,)).fetchall())
    below_count = project_count < current_app.config["PRIVATE_PROJECT_COUNT_LIMIT"]
    return below_count


def add_user_to_project(id_user: str,
                        id_project: str,
                        db,
                        role: str = "member",
                        permissions: dict = None):
    """
    Adds a user to a project. If permissions is "None", default is to add as data-read only.
    Parameters
    ----------
    id_user : str
        ID of the user to add.
    id_project : str
        ID of the project to which the user should be added.
    db : sqlite3.Connection
        Connection to the SQL database.
    role : str
        Name of the user's role.
    permissions : dict
        Optional. Dictionary keyed by permission name : value. See the SQL schema for possible keys.

    Returns
    -------
    None
    """
    if permissions is None:
        permissions = {"data_read": 1}
    elif 'data_read' not in permissions.keys() or permissions['data_read'] != 1:
        # Raise an error; caller is doing something weird.
        raise ValueError("The data_read property must be included and set to 1 for users of a project.")

    # Build command string
    cmd = f"INSERT INTO projects_roles (id_user, id_project, role"
    values = 'VALUES(?,?,?'
    for perm_name in permissions.keys():
        cmd += ',' + str(perm_name)
        values += ',?'
    cmd += ') '
    values += ')'
    cmd += values
    db.execute(cmd, (id_user, id_project, role) + tuple(permissions.values()))

    # Update number of users in project
    project_members = db.execute("SELECT role FROM projects_roles WHERE id_project = ?", (id_project,)).fetchall()
    db.execute("UPDATE projects SET num_members = ? WHERE id_project = ?", (len(project_members), id_project))
    db.commit()
    return


@bp.route('/results/<userid>/<filename>', methods=['GET'])
@login_required
def get_file(userid, filename):
    return send_from_directory(current_app.config['RESULTS_DIRECTORY'].format(userid=userid), filename)


# We are mixing several levels of abstraction; the programmer is being RUDE.
@bp.route('/project/<id_project>', methods=['GET', 'POST'])
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
            return url_for('pscs.project', id_project=id_project)
        elif 'delete' in request.json:  # is delete
            # Check that user is allowed
            has_perm = check_user_permission(permission_name='project_management',
                                             permission_value=1,
                                             id_project=id_project)
            if not has_perm:
                flash("You do not have permission to delete this project.")
                return url_for('pscs.project', id_project=id_project)
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
                    return url_for('pscs.project', id_project=id_project)
                add_user_to_project(id_user=id_user["id_user"], id_project=id_project, db=db,
                                    role='member', permissions={'data_read': 1})
            else:
                flash("You do not have permission to add users to the project.")
            return url_for("pscs.project", id_project=id_project)
        elif "deleteData" in request.json:
            id_data = request.json["deleteData"]
            delete_data(id_data)
            return url_for("pscs.project", id_project=id_project)
    return redirect(url_for('pscs.index'))


@bp.route('/project/<id_project>/public', methods=['GET', 'POST'])
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
    project_summary = db.execute("SELECT id_project, name_project, description "
                                 "FROM projects "
                                 "WHERE id_project = ?", (id_project,)).fetchone()
    # get results
    results_rows = db.execute(
        "SELECT file_path, title, description, id_analysis FROM results WHERE id_project = ?",
        (id_project,)).fetchall()
    results_files = []
    for r in results_rows:
        rr = dict(r)
        rr["file_name"] = os.path.basename(r["file_path"])
        results_files.append(rr)

    return render_template("pscs/project_public.html",
                           project_name=project_summary["name_project"],
                           summary=project_summary,
                           results=results_files)


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
        return render_template("pscs/project.html",
                               project_name=project_name,
                               analyses=analyses,
                               files=files,
                               analysis_nodes=analysis_nodes,
                               project_data_summary=project_data_summary,
                               results=results_files,
                               user_list=user_list,
                               summary=summary)


@bp.route('/project/<id_project>/publish', methods=["GET", "POST"])
@login_required
def project_publish(id_project):
    if request.method == "GET":
        has_perm = check_user_permission(permission_name='project_management',
                                         permission_value=1,
                                         id_project=id_project)
        if not has_perm:
            flash("You do not have permission to do this; contact the project owner.")
            return redirect(url_for("pscs.project", id_project=id_project))
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
        return render_template("pscs/project_publish.html",
                               proj=project_info,
                               authors=authors,
                               missing_name=missing_name,
                               project_data=project_data,
                               analyses=analyses,
                               results=result,
                               unrun_analyses=unrun_analyses)
    elif request.method == "POST":
        # Publish project
        # check perm
        has_perm = check_user_permission(permission_name="project_management",
                                         permission_value=1,
                                         id_project=id_project)
        if not has_perm:
            # Users will only end up here if their permissions have been changed while trying to publish a project
            # (or they're being rude)
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
            return redirect(url_for("pscs.project_publish", id_project=id_project))

        # Checks have been passed; user is allowed and has entered the correct confirmation.
        # Validate the user input.
        publication_info = request.json
        # Now need to check all data
        # Four things are sent: publication_type, authorlist, analyses, data
        # Validate publication type
        valid_pubtype = validate_publication_type(publication_info["publication_type"])
        valid_authors = validate_authorlist(publication_info["authorlist"], id_project=id_project)
        valid_analyses = validate_analyses(publication_info["analyses"], id_project=id_project)
        valid_data = validate_data(publication_info["data"], id_project=id_project)

        # Get author names
        if valid_pubtype and valid_authors and valid_analyses and valid_data:
            # Fetch results associated with analyses
            publication_info["results"] = {}
            for an_id in publication_info["analyses"]:
                publication_info["results"][an_id] = db.execute("SELECT id_result "
                                                                "FROM results "
                                                                "WHERE id_analysis = ? AND id_project = ?", (an_id, id_project))
            if publication_info["publication_type"] == "peer":
                publish_peer_review(id_project, project_summary=publication_info)
                flash(f"Project {name_project} set for peer review.")
            return_url = url_for('pscs.index')
            response_json = {"url": return_url, "publish_success": 1, "success_message": ""}
            return jsonify(response_json)
    return redirect(url_for("pscs.index"))


def validate_publication_type(pubtype: str) -> bool:
    return pubtype in ["peer", "public"]


def validate_authorlist(authorlist: list,
                        id_project: str) -> bool:
    """
    Checks whether each of the user IDs are associated with the project.
    Parameters
    ----------
    authorlist : list
        List of user IDs related to the project.
    id_project : str
        ID of the project that is being published.

    Returns
    -------
    bool
        Whether the list of authors is valid.
    """
    authorset = set(authorlist)
    # Get list of all authors associated with project; compare overlap
    db = get_db()
    authors = db.execute("SELECT id_user "
                         "FROM projects_roles "
                         "WHERE id_project = ?", (id_project,)).fetchall()
    author_ids = set([a["id_user"] for a in authors])
    return authorset == authorset.intersection(author_ids)


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
    missing_name = False
    for a in authors:
        if a["name"] is None or len(a["name"]) == 0:
            missing_name = True
            break
    return authors, missing_name


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


def publish_peer_review(id_project: str,
                        project_summary: dict):
    """
    Sets a project to peer review. Selected data, analyses, and results are placed behind a password prompt that does
    not require login.
    Parameters
    ----------
    id_project : str
        ID of the project to be published.
    project_summary : dict
        Dict containing what from the project should be made available.
    Returns
    -------
    None
    """
    # Need to do the following:
    db = get_db()
    # Mark project as under peer review
    db.execute("UPDATE projects "
               "SET is_peer_review = 1 "
               "WHERE id_project = ?", (id_project,))

    # Mark data to be under peer review.
    for id_data in project_summary["data"]:
        db.execute("UPDATE data "
                   "SET is_peer_review = 1 "
                   "WHERE id_data = ?", (id_data,))

    # Mark analysis to be under peer review
    for id_analysis in project_summary["analyses"]:
        db.execute("UPDATE analysis "
                   "SET is_peer_review = 1 "
                   "WHERE id_analysis = ?", (id_analysis,))

    # Mark results for each analysis to be under peer review.
    for id_an, results_list in project_summary["results"].items():
        for id_result in results_list:
            db.execute("UPDATE results "
                       "SET is_peer_review = 1 "
                       "WHERE id_result = ?", (id_result))

    # Store authors in order
    author_addresses = []
    for position, id_user in enumerate(project_summary["authorlist"]):
        db.execute("INSERT INTO publication_authors "
                   "(id_project, id_user, author_position) "
                   "VALUES (?,?,?)", (id_project, id_user, position))
        # Get author emails for later
        author_addresses.append(db.execute("SELECT email "
                                           "FROM users_auth "
                                           "WHERE id_user = ?", (id_user,)).fetchone()["email"])

    peer_review_password = generate_password()
    # Enter password into db
    pass_hash = generate_password_hash(peer_review_password)
    # We're checking the table twice; it's expected to be a small table, so we're going to rationalize our bad habits
    previous = db.execute("SELECT id_project "
                          "FROM projects_peer_review "
                          "WHERE id_project=?", (id_project,)).fetchone()
    if previous is None:
        db.execute("INSERT INTO projects_peer_review "
                   "(id_project, peer_password) "
                   "VALUES(?,?)", (id_project, pass_hash))
    else:
        # Overwrite previous record
        db.execute("UPDATE projects_peer_review "
                   "SET peer_password=? "
                   "WHERE id_project=?", (pass_hash, id_project))
    # Password is destined for HTML email; need to escape characters
    peer_review_password_html = werkzeug.utils.escape(peer_review_password)
    # Prepare email to notify user
    from pscs.templates.misc import peer_review_password_email
    project_info = db.execute("SELECT name_project, description "
                              "FROM projects "
                              "WHERE id_project = ?", (id_project,)).fetchone()
    project_url = "https://" + current_app.config["CURRENT_URL"] +  url_for("pscs.public_project", id_project=id_project)
    project_url = f"<a href='{project_url}'>{project_url}</a>"  # to make email content clickable
    peer_email = peer_review_password_email.format(title=project_info["name_project"],
                                                   description=project_info["description"],
                                                   project_url=project_url,
                                                   peer_review_password=peer_review_password_html)
    send_email(author_addresses, subject="PSCS Project Set for Peer Review", body=peer_email)
    db.commit()
    return


def generate_password(length: int = 32) -> str:
    """Generates a password"""
    charset = string.digits + string.ascii_letters + string.punctuation
    return ''.join(random.SystemRandom().choice(charset) for _ in range(length))


def check_user_permission(permission_name: str,
                          permission_value: int,
                          id_project: str) -> bool:
    """
    Checks that the user has the appropriate permission for the specified project.
    Parameters
    ----------
    permission_name : str
        Name of the permission to check.
    permission_value
        Value that the permission should have to be accepted.
    id_project : str
        ID of the project to check.

    Returns
    -------
    bool
        Whether the current user has permission.
    """
    db = get_db()
    role_info = db.execute(f'SELECT {permission_name} FROM projects_roles WHERE id_user = ? and id_project = ?',
                           (g.user['id_user'], id_project)).fetchone()
    if role_info is None:
        return False
    return role_info[permission_name] == permission_value


@bp.route('/run_analysis', methods=['POST'])
@login_required
def run_analysis():
    if request.method == 'POST':
        pipeline_specs = request.json
        id_project = session['CURRENT_PROJECT']
        id_analysis = pipeline_specs['id_analysis']
        # Get analysis json
        db = get_db()

        # Check whether user has already submitted a job
        id_user = g.user['id_user']
        file_ids = pipeline_specs['file_paths']

        submit_ok, reason = can_user_submit(id_user, id_project, id_analysis, file_ids, return_reason=True)
        if not submit_ok:
            return_url = url_for('pscs.project', id_project=session['CURRENT_PROJECT'])
            return jsonify({"url": return_url, "submit_status": reason, "submit_success": 0})

        # TODO: need to confirm that current user can access this analysis
        pipeline_json = db.execute('SELECT node_file FROM analysis WHERE id_analysis = ?',
                                   (pipeline_specs['id_analysis'],)).fetchall()[0]['node_file']

        # This next section gets the relevant paths for input and output.
        # Instead, it should create the output paths (as needed), gather the input files, and trigger for them to be
        # sent out to OSG.
        output_dir = current_app.config['RESULTS_DIRECTORY'].format(id_project=session['CURRENT_PROJECT'],
                                                            id_analysis=pipeline_specs['id_analysis'])
        pathlib.Path(output_dir).mkdir(exist_ok=True)
        file_ids = pipeline_specs['file_paths']
        file_info = dict()
        for node_id, file_id in file_ids.items():
            buf = db.execute('SELECT file_path FROM data WHERE id_data = ?', (file_id,)).fetchall()[0]
            file_info[node_id] = dict()
            file_info[node_id]["file_path"] = buf["file_path"]
            file_info[node_id]["id"] = file_id

        # Dispatch to OSP
        dispatch(pipeline_json=pipeline_json,
                 id_user=g.user['id_user'],
                 file_info=file_info,
                 id_project=id_project,
                 id_analysis=id_analysis,
                 resource='osp')
        response_json = {"url": url_for('pscs.project', id_project=session['CURRENT_PROJECT']),
                         "redirect": True}
        return jsonify(response_json)


def delete_data(id_data):
    """
    Deletes the data specified by the POST request. Verifies that user has permission to do so.
    Returns
    -------
    None
    """
    # Get related project id
    db = get_db()
    data_row = db.execute('SELECT * FROM data WHERE id_data = ?', (id_data,)).fetchone()
    data_write = check_user_permission(permission_name="data_write",
                                       permission_value=1,
                                       id_project=data_row["id_project"])
    if not data_write:
        flash("You do not have permission to delete data. Contact your project's manager to remove data.")
    elif data_row['is_published']:
        flash("Data has been published and can't be deleted.")
    elif data_write and data_row['is_published'] == 0:  # assertion of conditions instead of 'else'
        # Permissions check out; stage data for deletion
        db.execute('DELETE FROM data WHERE id_data = ?', (id_data,))
        db.commit()
        flash("Data deleted.")
    return


@bp.route('/pipeline', methods=['GET', 'POST'])
def pipeline_designer():
    if request.method == 'GET':
        static_path = current_app.config["STATIC_DIRECTORY"]
        f = open(os.path.join(static_path, "node_data.json"), 'r')
        j = json.load(f)
        f.close()
        modules = get_unique_values_for_key(j, 'module')
        node_json = convert_dict_to_list(j)
        proj_dests, user_dests = get_save_destinations()

        # Fetch analysis that user is allowed to read
        db = get_db()
        analyses = db.execute("SELECT analysis.id_analysis, analysis.analysis_name "
                              "FROM analysis INNER JOIN projects_roles "
                              "ON analysis.id_project=projects_roles.id_project "
                              "WHERE projects_roles.id_user = ? AND projects_roles.analysis_read = 1",
                              (g.user['id_user'],)).fetchall()
        if "CURRENT_PROJECT" not in session.keys():
            return render_template("pscs/pipeline.html", node_json=node_json, modules=modules, proj_dests=proj_dests,
                                   user_dests=user_dests, analyses=analyses)
        else:
            return render_template("pscs/pipeline.html", node_json=node_json, modules=modules, proj_dests=proj_dests,
                                   user_dests=user_dests, analyses=analyses, current_project=url_for("pscs.project", id_project=session["CURRENT_PROJECT"]))
    elif request.method == 'POST':
        if "loadAnalysis" in request.json:
            # Check user perm
            # Get id of related project
            db = get_db()
            id_analysis = request.json["loadAnalysis"]
            id_project = db.execute("SELECT id_project "
                                    "FROM analysis WHERE id_analysis = ?", (id_analysis,)).fetchone()["id_project"]
            has_perm = check_user_permission(permission_name="analysis_read",
                                             permission_value=1,
                                             id_project=id_project)
            if not has_perm:
                return {"": ""}
            elif has_perm:
                # has permission to read; go get analysis file and return JSON
                node_file = db.execute("SELECT node_file FROM analysis WHERE id_analysis = ?", (id_analysis,)).fetchone()['node_file']
                f = open(node_file, 'r')
                node_data = json.load(f)
                f.close()
                return node_data
            return {"": ""}
        pipeline_summary = request.json
        is_dest_project = pipeline_summary['isDestProject']
        id_project = pipeline_summary['saveId']
        # Go through nodes to find which ones have path keyword as arguments; this identifies inputs
        input_nodes = {}
        for n in pipeline_summary['nodes']:
            params = n['paramsValues']
            if PATH_KEYWORD in params.keys():
                input_nodes[n['nodeId']] = n['labelText']  # labelText is what is displayed to the user
        # TODO: saveID is from the user; need to validate that the user has access to it
        # Create new pipeline ID
        db = get_db()
        id_analysis = get_unique_value_for_field(db, "id_analysis", "analysis")
        pipeline_id = id_analysis

        if is_dest_project:
            pipeline_dir = current_app.config['PROJECTS_DIRECTORY'].format(id_project=id_project)
        # else:
            # pipeline_dir = current_app.config['UPLOAD_FOLDER'].format(userid=g['id_user'])

        output_name = pipeline_id + '.json'
        pipeline_file = os.path.join(pipeline_dir, output_name)
        f = open(pipeline_file, 'w')
        json.dump(pipeline_summary, f, indent=2)
        f.close()

        pipeline_hash = calc_hash_of_file(pipeline_file)
        pipe_name = secure_filename(pipeline_summary['name'])
        # send to database
        if is_dest_project:
            db.execute(
                'INSERT INTO analysis '
                '(id_analysis, id_project, analysis_name, node_file, parameter_file, analysis_hash)'
                ' VALUES (?,?,?,?,?,?)',
                (id_analysis, id_project, pipe_name, pipeline_file, pipeline_file, pipeline_hash))
            for inp_id, inp_name in input_nodes.items():
                # get uuid
                id_input = get_unique_value_for_field(db, 'id_input', 'analysis_inputs')
                db.execute(
                    'INSERT INTO analysis_inputs (id_input, id_analysis, node_id, node_name)'
                    ' VALUES (?,?,?,?)', (id_input, id_analysis, inp_id, inp_name)
                )
            db.commit()
        return render_template("pscs/pipeline.html")


@bp.route("/about", methods=["GET"])
def about():
    if request.method == "GET":
        return render_template("about/about.html")


def get_save_destinations():
    # Get all projects that this user belongs to
    db = get_db()
    projs = db.execute(
        'SELECT name_project, id_project FROM projects WHERE id_user = ?',
        (g.user['id_user'],)).fetchall()
    proj_dests = {}
    for p in projs:
        proj_dests[p['name_project']] = p['id_project']
    user_dests = {'user': g.user['name_user'], 'id': g.user['id_user']}
    # user_dests = {g.user['name_user']: g.user['id_user']}
    # return [p['name_project'] for p in projs]
    return proj_dests, user_dests


def convert_dict_to_list(d: dict) -> list:
    """
    Converts a dict of dicts into a list of dicts with the initial key converted to "name".
    This function is intended to be used to pass JSON contents to Flask.
    Parameters
    ----------
    d : dict
        Dict of dicts to convert.
    Returns
    -------
    list
        List of dicts.
    """
    dict_list = []
    for k, v in d.items():
        v['name'] = k
        dict_list.append(v)
    return dict_list


def get_unique_values_for_key(d: dict, key) -> list:
    """
    Gets the unique values for a given key.
    Parameters
    ----------
    d : dict
        Dict to check
    key
        Key to use for value lookups.

    Returns
    -------
    list
        List of unique values
    """
    unique_values = set()
    for k, v in d.items():
        unique_values.add(v[key])
    return sorted(unique_values)


def calc_hash_of_file(file: str) -> str:
    """
    Calculates the SHA256 hash of a file by loading 1kB at a time.
    Parameters
    ----------
    file : str
        Path to the file for which to compute the hash.

    Returns
    -------
    str
        SHA256 hash of the file.
    """
    sha = hashlib.sha256()
    f = open(file, 'rb')
    dat = f.read(1024)
    while dat:  # iterate until end of file
        sha.update(dat)
        dat = f.read(1024)
    return sha.hexdigest()


@bp.route('/projects/<id_project>/results/<id_analysis>/<path:filename>', methods=['GET'])
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
