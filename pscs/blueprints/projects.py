import sqlite3

from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app, send_from_directory,
    jsonify, send_file
)
from collections import defaultdict as dd
import datetime
import werkzeug
import secrets
from werkzeug.utils import secure_filename, escape
import email_validator
from email_validator import EmailNotValidError
from pscs.auth import login_required, is_logged_in
from pscs.db import get_db, check_user_permission, check_analysis_published, get_default_analysis
from pscs.pscs import delete_data, get_unique_value_for_field, add_user_to_project
from pscs.transfers.fetching import read_logs
from werkzeug.security import check_password_hash, generate_password_hash
from pscs.messaging.mail import send_email
import pscs
from pscs.pscs import load_analysis_from_id
import os
from os.path import join
import shutil
import string
import random
from itsdangerous.url_safe import URLSafeTimedSerializer
from typing import Collection, Optional
import zipfile
import tempfile
from pscs.extensions.limiter import limiter
import requests
import time
from pybtex.database.input import bibtex
import pickle as pkl


bp = Blueprint("projects", __name__, url_prefix="/project")


@bp.route('/<id_project>', methods=['GET', 'POST'])
@login_required
def project(id_project):
    if request.method == 'GET':
        return display_private_project(id_project)
    elif request.method == 'POST':
        # Check for rename or delete
        if 'newName' in request.json:  # is rename
            new_project_name = escape(request.json['newName'])
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
        elif 'inviteUser' in request.json:
            # check that current user is allowed to add users
            has_perm = check_user_permission(permission_name='project_management',
                                             permission_value=1,
                                             id_project=id_project)
            if has_perm:
                # add user
                db = get_db()
                user_to_invite = escape(request.json['inviteUser'])
                # Get user id from email or username
                try:
                    email_validator.validate_email(user_to_invite)
                    id_user = db.execute("SELECT id_user FROM users_auth WHERE email = ?", (user_to_invite,)).fetchone()
                except EmailNotValidError:
                    # not valid email; done by username
                    id_user = db.execute("SELECT id_user FROM users_auth WHERE name_user = ?", (user_to_invite,)).fetchone()

                if id_user is None:
                    # user doesn't exist
                    flash(f"User {user_to_invite} not found.")
                    return url_for('projects.project', id_project=id_project)
                invite_user_to_project(id_invitee=id_user["id_user"], id_inviter=session["id_user"], id_project=id_project)
                flash("User invited.")
                return url_for("projects.project", id_project=id_project)
            else:
                flash("You do not have permission to add users to the project.")
            return url_for("projects.project", id_project=id_project)
        elif "deleteData" in request.json:
            id_data = request.json["deleteData"]
            has_perm = check_user_permission(permission_name="data_write",
                                             permission_value=1,
                                             id_project=id_project)
            if has_perm:
                delete_data(id_data)
            else:
                flash("You do not have permission to delete data from the project.")
            return url_for("projects.project", id_project=id_project)
        elif "description" in request.json:
            has_perm = check_user_permission(permission_name="project_management",
                                             permission_value=1,
                                             id_project=id_project)
            if not has_perm:
                flash("You do not have permission to change the project description.")
                return url_for("projects.project", id_project=id_project)
            # First escape description
            new_description = escape(request.json["description"])
            if len(new_description) > current_app.config["MAX_DESCRIPTION_LENGTH"]:
                    new_description = new_description[:current_app.config["MAX_DESCRIPTION_LENGTH"]] + " [...]"
                    flash(f"The project description exceeds the maximum length ({current_app.config['MAX_DESCRIPTION_LENGTH']}). It has been shortened.")
            db = get_db()
            db.execute("UPDATE projects SET description = ? WHERE id_project = ?", (new_description, id_project))
            db.commit()
            return url_for("projects.project", id_project=id_project)

    return redirect(url_for('pscs.index'))


@bp.route("/manage_invitation", methods=["POST"])
def manage_invitation():
    if request.method == "POST":
        data = request.json
        if "invite" == data["action"]:
            if check_user_permission("project_management", 1, data["id_project"]):
                invite_user_to_project(id_invitee=data["id_invitee"],
                                       id_inviter=session["id_user"],
                                       id_project=data["id_project"])
                return jsonify({"done": 1})
        elif "rescind" == data["action"]:
            # Confirm that user either sent the invitation (somehow) or has management permission
            invitation = get_invitation(data["id_invitation"])
            if invitation["id_inviter"] == session["id_user"] or \
                    check_user_permission("project_management", 1, invitation["id_project"]):
                rescind_invitation(data["id_invitation"])
                flash("Invitation rescinded.")
                reply_json = {"done": 1, "url": url_for("pscs.projects_summary")}
                return jsonify(reply_json)
            flash("You do not have the permission to do this.")
            return redirect(url_for("pscs.index"))
        elif "accept" == data["action"]:
            # Confirm invitation is real
            db = get_db()
            invitation = get_invitation(data["id_invitation"])
            if invitation["id_invitee"] == session["id_user"]:
                # invitation exists, etc.
                add_user_to_project(id_user=session["id_user"],
                                    id_project=invitation["id_project"],
                                    db=db,
                                    role="member",
                                    permissions={"data_read": 1,
                                                 "data_write": 1,
                                                 "analysis_read": 1,
                                                 "analysis_write": 1,
                                                 "analysis_execute": 1,
                                                 "project_management": 0})
                rescind_invitation(data["id_invitation"])
                return jsonify({"done": 1})
        elif "reject" == data["action"]:
            invitation = get_invitation(data["id_invitation"])
            if invitation["id_invitee"] == session["id_user"]:
                # invitation is for the current user
                rescind_invitation(data["id_invitation"])
                return jsonify({"done": 1, "url": url_for("pscs.projects_summary")})
    return redirect(url_for("pscs.index"))


@bp.route("/<id_project>/upload", methods=["GET"])
@login_required
def upload(id_project):
    if request.method == "GET":
        return render_template("pscs/uploader/landing.html")


def get_filestorage_size(f) -> int:
    file_length = f.seek(0, os.SEEK_END)
    f.seek(0)
    return file_length


def get_invitation(id_invitation: str):
    """Fetches inviter, invitee, and project ids connected to the invitation"""
    db = get_db()
    return db.execute("SELECT id_invitee, id_inviter, id_project "
                      "FROM project_invitations "
                      "WHERE id_invitation = ?", (id_invitation,)).fetchone()


def invite_user_to_project(id_invitee: str,
                           id_inviter: str,
                           id_project: str):
    """
    Sends invitation to the specified user to join the project.
    Parameters
    ----------
    id_invitee : str
        ID of the user to invite.
    id_inviter : str
        ID of the user doing the inviting.
    id_project : str
        ID of the project.

    Returns
    -------
    None
    """
    db = get_db()
    id_invitation = get_unique_value_for_field(db, field="id_invitation", table="project_invitations")
    db.execute("INSERT INTO project_invitations (id_invitation, id_invitee, id_inviter, id_project) VALUES(?,?,?,?)",
               (id_invitation, id_invitee, id_inviter, id_project))
    db.commit()
    return


def rescind_invitation(id_invitation: str):
    """
    Rescinds a previously-sent invitation.
    Parameters
    ----------
    id_invitation : str
        ID of the invitation to rescind.

    Returns
    -------
    None
    """
    db = get_db()
    db.execute("DELETE FROM project_invitations "
               "WHERE id_invitation = ?", (id_invitation,))
    db.commit()
    return


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
            return redirect(url_for("projects.public_project", id_project=id_project))
    return


@bp.route('/<id_project>/public/load_analysis', methods=["POST"])
def load_analysis(id_project):
    if request.method == "POST":
        if "loadAnalysis" in request.json:
            id_analysis = request.json["loadAnalysis"]
            # Verify that the project is indeed public
            db = get_db()
            status = _get_project_publication_status(db, id_project)
            has_peer_review_password = status == "peer review" and "project_review" in session.keys() and id_project in session["project_review"]
            if status == "public" or has_peer_review_password:
                return load_analysis_from_id(id_analysis)
            elif check_user_permission("analysis_read", 1, id_project):
                # has permission to read; go get analysis file and return JSON
                return load_analysis_from_id(id_analysis)
            return {"": ""}


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
                          "FROM analysis AS A "
                          "WHERE A.id_project = ? AND (is_peer_review=1 OR is_published=1)",
                          (id_project,)).fetchall()
    project_summary["analyses"] = [dict(an) for an in analyses]
    public_analysis_id = set()
    for an in project_summary["analyses"]:
        public_analysis_id.add(an["id_analysis"])
    results_rows = db.execute("SELECT file_path, file_name, title, description, id_analysis "
                              "FROM results "
                              "WHERE id_project = ? AND (is_peer_review=1 OR is_published=1)", (id_project,)).fetchall()
    # Single DB fetch op means we have to sort them by analysis now
    project_summary["results"] = dd(list)
    for r in results_rows:
        an_id = r["id_analysis"]
        project_summary["results"][an_id].append(dict(r))
        project_summary["results"][an_id][-1]["file_name"] = os.path.basename(project_summary["results"][an_id][-1]["file_path"])
        project_summary["results"][an_id][-1]["id_job"] = os.path.split(project_summary["results"][an_id][-1]["file_path"])[-1]

    # For logged-in users, make it easy to import pipelines
    # Get user's projects
    user_projects = None
    if "id_user" in session and session["id_user"] is not None:
        user_projects = db.execute("SELECT P.id_project, P.name_project "
                                   "FROM projects AS P INNER JOIN "
                                   "projects_roles AS PR ON P.id_project = PR.id_project WHERE "
                                   "PR.id_user = ? AND PR.analysis_write = 1", (session["id_user"],)).fetchall()

    # For projects that are under peer-review, or for public projects where the user is logged in, make it possible
    # to download the published data.
    project_status = _get_project_publication_status(db, id_project=id_project)
    project_data_list = None
    analysis_list = None
    if (project_status == "peer_review") or (project_status == "public" and "id_user" in session and session["id_user"] is not None):
        # Get list of data associated with the project
        project_data = db.execute("SELECT D.id_data, D.file_name "
                                  "FROM data AS D INNER JOIN projects AS P "
                                  "ON D.id_project = P.id_project "
                                  "WHERE D.id_project = ? AND (D.is_published = 1 OR D.is_peer_review = 1)", (id_project,)).fetchall()
        project_data_list = []
        for pd in project_data:
            project_data_list.append({"id_data": pd["id_data"], "file_name": pd["file_name"]})

        # Get analysis for downloading results
        analysis_data = db.execute("SELECT id_analysis, analysis_name "
                                   "FROM analysis WHERE id_project=? AND (is_published = 1 OR is_peer_review = 1)", (id_project,)).fetchall()
        analysis_list = []
        for an in analysis_data:
            analysis_list.append({"id_analysis": an["id_analysis"], "analysis_name": an["analysis_name"]})

    return render_template("pscs/project_public.html",
                           project_summary=project_summary,
                           user_projects=user_projects,
                           project_data=project_data_list,
                           project_analysis=analysis_list)


@bp.route('/<id_project>/results/<id_analysis>/<id_job>/<path:filename>', methods=['GET'])
def results(filename, id_project, id_job, id_analysis):
    # If logged in and has permission, don't need to check public
    if is_logged_in() and check_user_permission("data_read", 1, id_project):
        return private_results(filename, id_project, id_analysis)
    db = get_db()
    id_result = filename.split(os.path.extsep)[0]
    public_status = db.execute("SELECT P.status FROM publications AS P INNER JOIN publications_results AS PR ON "
                               "P.id_publication = PR.id_publication WHERE PR.id_result = ?", (id_result,)).fetchone()
    if public_status is None:
        return  # Problem
    if public_status["status"] == "public":
        return public_results(filename, id_project, id_analysis)
    elif public_status["status"] == "peer review":
        return review_results(filename, id_project, id_analysis)


# @login_required
def private_results(filename, id_project, id_analysis):
    if is_logged_in():
        has_perm = check_user_permission("data_read", 1, id_project)
        if not has_perm:
            return
        id_job = filename.split(os.path.extsep)[0]
        res_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis, id_job=id_job)
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
        default_analyses = get_default_analysis()
        for an in default_analyses:
            analyses[an["id_analysis"]] = an["analysis_name"]
            input_nodes = db.execute("SELECT node_id, node_name FROM analysis_inputs WHERE id_analysis = ?", (an["id_analysis"],)).fetchall()
            project_inps = {}
            for inp in input_nodes:
                project_inps[inp["node_id"]] = inp["node_name"]
            analysis_nodes[an["id_analysis"]] = project_inps

        # Get files associated with project
        project_data = db.execute('SELECT id_data, file_name FROM data WHERE id_project = ?',
                                  (id_project,)).fetchall()
        files = {}
        for project_file in project_data:
            files[project_file['id_data']] = project_file["file_name"]
        project_data_summary = db.execute(
            'SELECT id_data, file_path, data_type, file_hash, data_uploaded_time FROM data WHERE id_project = ?',
            (id_project,)).fetchall()
        # Convert returned data to dicts; remove leading path
        new_summary = []
        for dat in project_data_summary:
            ddat = dict(dat)
            ddat["file_path"] = os.path.basename(ddat["file_path"])
            new_summary.append(ddat)
        project_data_summary = new_summary

        # Get results
        results_rows = db.execute("SELECT file_path, title, description, id_analysis FROM results WHERE id_project = ?",
                                  (id_project,)).fetchall()
        results_files = []
        for r in results_rows:
            rr = dict(r)
            rr["file_path"] = rr["file_path"][len(current_app.config["INSTANCE_PATH"]):]
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
        public_status = _get_project_publication_status(db, id_project)

        project_summary = get_project_summary(db, id_project)
        project_summary["id_project"] = id_project
        return render_template("pscs/project.html",
                               project_name=project_name,
                               analyses=analyses,
                               files=files,
                               analysis_nodes=analysis_nodes,
                               project_data_summary=project_data_summary,
                               results=results_files,
                               user_list=user_list,
                               project_status=public_status,
                               summary=project_summary)
    else:
        flash("The requested project is not available.")
        return redirect(url_for("pscs.index"))


def _get_publication_project(id_publication):
    db = get_db()
    return db.execute("SELECT id_project FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()["id_project"]


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
        public_status = _get_project_publication_status(db, id_project)
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
        user_submmited_name = escape(request.json["confirmation"])
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
        if all_valid:
            # Fetch results associated with analyses
            publication_info["results"] = {}
            for an_id in publication_info["analyses"]:
                buf = db.execute("SELECT id_result "
                                 "FROM results "
                                 "WHERE id_analysis = ? AND id_project = ?", (an_id, id_project)).fetchall()
                publication_info["results"][an_id] = [b["id_result"] for b in buf]
            db = get_db()
            has_external_authors = _contains_external_authors(publication_info["authorlist"])
            # Three possibilities: public, peer review, update
            if not has_external_authors:
                publication_info["hold"] = False
                id_publication = prepare_publication(id_project=id_project, publication_info=publication_info, contact_user=g.user['id_user'], commit=True)
                return_url = url_for("publications.request_publication", id_publication=id_publication)
                response_json = {"url": return_url, "publish_success": 1, "success_message": ""}
                if publication_info["publication_type"] == "public":
                    notify_all_authors_public(id_publication)
                    flash(f"The project has been set to \"{publication_info['publication_type']}\"")
                    return jsonify(response_json)
                elif publication_info["publication_type"] == "peer review":
                    # Only send password to publishing author.
                    current_info = _get_publication_info(db, id_project=id_project)[-1]  # only get latest
                    peer_review_password = generate_peer_review_password(id_publication=id_publication)
                    notify_user_peer_review_password(id_user=session["id_user"],
                                                     publication_info=current_info,
                                                     peer_review_password=peer_review_password)
                    notify_authors_peer_review(publication_info=current_info,
                                               publishing_id_user=session["id_user"])
                    flash(f"The project has been set to \"{publication_info['publication_type']}\"")
                    return jsonify(response_json)
                elif publication_info["publication_type"] == "update":
                    return jsonify(response_json)
            else:
                # Project has external authors; must be placed on hold until they confirm
                publication_info["hold"] = True
                id_publication = prepare_publication(id_project=id_project, publication_info=publication_info, contact_user=g.user["id_user"], commit=False)
                return_url = url_for("pscs.index")
                response_json = {"url": return_url, "publish_success": 1, "success_message": "External authors must approve their association before proceeding."}
                notify_external_authors(id_publication=id_publication)
                flash(f"External authors have been notified; they will need to complete the information request before proceeding.")
                db.commit()
                return jsonify(response_json)

        else:
            # These only happen if the project changes while the user is trying to publish it, or sneaky user.
            flash_msg = []
            if not valid_pubtype:
                flash_msg.append(f"Project can't be published as {publication_info['publication_type']}")
            if not valid_authors:
                flash_msg.append("Some of the authors are PSCS users that are not associated with the project.")
            if not valid_analyses:
                flash_msg.append("Some of the selected analyses can't be published. Ensure that they have been run first.")
            if not valid_data:
                flash_msg.append("Some of the selected data is not associated with the project.")
            flash_msg = " - ".join(flash_msg)
            flash(flash_msg)
            return_url = url_for('pscs.index')
            response_json = {"url": return_url, "publish_success": 0, "success_message": flash_msg}
            return jsonify(response_json)
    return redirect(url_for("pscs.index"))


def _get_project_publication_status(db, id_project) -> str:
    """Returns the publication status of the specified project. Returns either "public", "peer review", or "private"."""
    status = db.execute("SELECT status "
                      "FROM publications "
                      "WHERE id_project = ? "
                      "ORDER BY creation_time DESC", (id_project,)).fetchone()
    if status is None:
        return "private"
    else:
        return status["status"]


def validate_publication_type(pubtype: str, id_project: str) -> bool:
    """Returns whether the project can be set to the specified publication type."""
    db = get_db()
    public_status = _get_project_publication_status(db, id_project)
    valid_pubtypes = ["peer review", "public", "update"]
    if pubtype == public_status:
        return False  # not doing anything
    elif pubtype not in valid_pubtypes:
        return False  # bad string
    elif public_status == "private":
        return True  # if it's private, any publication is fine (this one shouldn't happen)
    elif public_status == "peer review" and pubtype == "public":
        return True  # if it's under peer review, only 'public' is valid
    elif public_status == "public" and pubtype == "update":
        return True  # currently public, but is getting a new version
    elif pubtype == "update" and public_status != "public":
        return False  # only updates to public projects allowed
    else:
        # There is the rare circumstance where projects that are published or set for peer review need to be private
        # Users can do this only within the window specified in the data use agreement.
        if pubtype == "private":
            time_published = _get_project_publication_time(db, id_project)
            now_utc = int(datetime.datetime.now().astimezone(datetime.timezone.utc).strftime("%s"))
            return now_utc - time_published <= current_app.config["PUBLIC_RECANT_TIMEOUT"]
    return False  # default


def _get_project_publication_time(db, id_project: str) -> float:
    """Returns the epoch time for when a project was published."""
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
    try:
        author_list.sort(key=lambda x: x["author_position"])
    except KeyError:
        pass
    missing_name = False
    for a in authors:
        if a["name"] is None or len(a["name"]) == 0:
            missing_name = True
            break
    return author_list, missing_name


def _get_data(id_project) -> list:
    """Returns the list of data that has been uploaded to the project."""
    db = get_db()
    return db.execute("SELECT id_data, file_name, data_uploaded_time "
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
                        publication_info: dict,
                        contact_user: str,
                        commit: bool = False) -> str:
    """
    Sets a project to peer review or public. Selected data, analyses, and results are placed behind a password prompt that does
    not require login.
    Parameters
    ----------
    id_project : str
        ID of the project to be published.
    publication_info : dict
        Dict containing what from the project should be made available.
    commit : bool
        Whether to commit changes to database.
    Returns
    -------
    None
    """
    db = get_db()
    # Get title, etc. from project being published
    project_info = db.execute("SELECT name_project, description FROM projects WHERE id_project = ?", (id_project,)).fetchone()
    if publication_info["publication_type"] == "update":
        status = "public"
    else:
        status = publication_info["publication_type"]
    # Get current status
    current_publication = _get_publication_info(db, id_project)  # get latest
    if current_publication is not None and len(current_publication) > 0 and current_publication[-1]["status"] == "peer review" and status == "public":
        # Project is being moved from peer review to public
        db.execute("UPDATE publications SET status = ? WHERE id_publication = ?", (status, current_publication[-1]["id_publication"]))
        if commit:
            db.commit()
        return current_publication[-1]["id_publication"]
    # Generate publication ID
    id_publication = get_unique_value_for_field(db, "id_publication", "publications")
    db.execute("INSERT INTO publications (id_publication, id_project, title, description, status, hold, id_user_contact) "
               "VALUES (?,?,?,?,?,?,?)", (id_publication, id_project, project_info["name_project"], project_info["description"], status, publication_info["hold"], contact_user))

    # Move data into publication
    for id_data in publication_info["data"]:
        db.execute("INSERT INTO publications_data (id_publication, id_data) VALUES (?,?)", (id_publication, id_data))
    # Move analysis into publication
    for id_analysis in publication_info["analyses"]:
        db.execute("INSERT INTO publications_analysis (id_publication, id_analysis) VALUES (?,?)", (id_publication, id_analysis))
    # Move results into publication
    for id_analysis, results_list in publication_info["results"].items():
        for id_result in results_list:
            db.execute("INSERT INTO publications_results (id_publication, id_result) VALUES (?,?)", (id_publication, id_result))
    pscs_authors, external_authors = _split_authlist(publication_info["authorlist"])
    # Set authors for publication
    for author in pscs_authors:
        db.execute("INSERT INTO publications_authors (id_publication, id_user, author_position, confirmed) VALUES (?,?,?,?)",
                   (id_publication, author["id"], author["position"], 1))
        # For each author, get affiliation
        author_affiliation = db.execute("SELECT affiliation, affiliation_order FROM users_affiliation WHERE id_user = ?", (author["id"],)).fetchall()
        for aff in author_affiliation:
            db.execute("INSERT INTO publications_authors_affiliation (id_publication, id_user, affiliation, affiliation_order) VALUES(?,?,?,?)",
                       (id_publication, author["id"], aff["affiliation"], aff["affiliation_order"]))

    for author in external_authors:
        db.execute("INSERT INTO publications_external_authors_info (id_publication, email, author_position) VALUES (?,?,?)",
                   (id_publication, author["email"], author["position"]))


    # Make shortid
    _ = make_shortid(db, id_publication)
    if commit:
        db.commit()
    return id_publication


def get_short_url(id_publication):
    db = get_db()
    id_pub_short = db.execute("SELECT id_publication_short FROM publications_shortid AS PS WHERE PS.id_publication = ?", (id_publication,)).fetchone()["id_publication_short"]
    return build_full_url(url_for("publications_short.short_publication", id_publication_short=id_pub_short))


def notify_user_peer_review_password(id_user: str,
                                    publication_info: dict,
                                    peer_review_password: str):
    """Notifies the specified user that a publication has been set to peer review and gives them the password for submission."""
    db = get_db()
    user_email = db.execute("SELECT email FROM users_auth WHERE id_user = ?", (id_user,)).fetchone()["email"]
    from pscs.templates.misc import peer_review_password_email
    publication_url = build_full_url(url_for("publications.request_publication", id_publication=publication_info["id_publication"]))
    publication_shorturl = get_short_url(publication_info["id_publication"])

    peer_review_password_html = werkzeug.utils.escape(peer_review_password)

    peer_email = peer_review_password_email.format(title=publication_info["title"],
                                                   description=publication_info["description"],
                                                   publication_url=publication_url,
                                                   publication_shorturl=publication_shorturl,
                                                   peer_review_password=peer_review_password_html)
    send_email([user_email], "PSCS Project Set for Peer Review", body=peer_email)
    return


def notify_authors_peer_review(publication_info: dict,
                               publishing_id_user: str):
    """NOtifies authors associated with a publication that it has been set for peer review. The publishing author is
    excluded since they should receive the notification email containing the password."""
    emails = set(get_publication_author_emails(publication_info["id_publication"]))
    db = get_db()
    publishing_user_info = db.execute("SELECT name, name_user, email FROM users_auth WHERE id_user = ?", (publishing_id_user,)).fetchone()
    emails = emails.difference(set([publishing_user_info['email']]))
    if len(emails) == 0:
        return  # if there is only one author
    publication_url = build_full_url(url_for("publications.request_publication", id_publication=publication_info["id_publication"]))
    publication_shorturl = get_short_url(publication_info["id_publication"])

    from pscs.templates.misc import peerreview_notification
    peer_email = peerreview_notification.format(title=publication_info["title"],
                                                description=publication_info["description"],
                                                publication_url=publication_url,
                                                publication_shorturl=publication_shorturl,
                                                publishing_user=publishing_user_info["name_user"],
                                                publishing_name=publishing_user_info["name"])
    send_email(emails, "PSCS Project Set for Peer Review", body=peer_email)
    return


def notify_all_authors_public(id_publication):
    """
    Publishes the specified project to be publicly viewable.
    Parameters
    ----------
    id_publication : str
        ID of the associated publication.
    Returns
    -------
    None
    """
    db = get_db()
    # Format notification email
    from pscs.templates.misc import publication_notification
    id_project = db.execute("SELECT id_project FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()["id_project"]
    publication_info = _get_publication_info(db, id_project)[-1]  # get last publication
    publication_url = build_full_url(url_for("publications.request_publication", id_publication=id_publication))
    publication_shorturl = get_short_url(id_publication=id_publication)
    notif_email = publication_notification.format(title=publication_info["title"],
                                                  publication_url=publication_url,
                                                  publication_shorturl=publication_shorturl)
    author_addresses = get_publication_author_emails(id_publication)
    send_email(author_addresses, subject="PSCS Project Publication", body=notif_email)
    db.commit()
    return


def generate_peer_review_password(id_publication: str):
    peer_review_password = generate_password(length=32)
    pass_hash = generate_password_hash(peer_review_password)
    db = get_db()
    db.execute("INSERT INTO publications_peer_review (id_publication, peer_password) VALUES(?,?)", (id_publication, pass_hash))
    db.commit()
    return peer_review_password

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
    project_url = current_app.config["CURRENT_URL"] + url_for("projects.public_project", id_project=id_project)
    project_url = f"<a href='{project_url}'>{project_url}</a>"  # to make email content clickable
    peer_email = peer_review_password_email.format(title=project_info["name_project"],
                                                   description=project_info["description"],
                                                   project_url=project_url,
                                                   peer_review_password=peer_review_password_html)
    author_addresses = get_publication_author_emails(id_project)
    send_email(author_addresses, subject="PSCS Project Set for Peer Review", body=peer_email)
    db.commit()
    return





def _remove_peer_review_password(id_project, commit=False):
    """Removes the project's password associated with peer review, if any."""
    db = get_db()
    db.execute("DELETE FROM projects_peer_review WHERE id_project = ?", (id_project,))
    if commit:
        db.commit()
    return


def generate_password(length: int = 32,
                      charset: str = None) -> str:
    """Generates a password of the specified length. Defaults"""
    if charset is None:
        charset = string.digits + string.ascii_letters
    return ''.join(secrets.choice(charset) for _ in range(length))


# NOTE: These are constructed as separate functions to avoid parameterizing the SQL string.
def _set_project_status(id_project, is_peer_review=None, is_published=None, project_version=0, commit=False):
    """Sets the publication status bits (is_peer_review, is_published) for the specified project. If None, their value
    is unchanged."""
    db = get_db()
    if is_peer_review is not None and is_published is not None:
        db.execute("UPDATE projects "
                   "SET is_peer_review = ?, is_published = ?, version = ? "
                   "WHERE id_project = ?", (is_peer_review, is_published, project_version, id_project))
    elif is_peer_review is not None:
        db.execute("UPDATE projects SET is_peer_review = ? WHERE id_project = ?", (is_peer_review, id_project))
    elif is_published is not None:
        db.execute("UPDATE projects SET is_published = ?, version = ? WHERE id_project = ?", (is_published, project_version, id_project))
    if commit:
        db.commit()
    return


def _set_data_status(id_data, is_peer_review=None, is_published=None, project_version=0, commit=False):
    """Sets the publication status bits (is_peer_review, is_published) for the specified datum. If None, their value is
    unchanged."""
    db = get_db()
    if is_peer_review is not None and is_published is not None:
        db.execute("UPDATE data "
                   "SET is_peer_review = ?, is_published = ?, project_version = ? "
                   "WHERE id_data = ?", (is_peer_review, is_published, project_version, id_data))
    elif is_peer_review is not None:
        db.execute("UPDATE projects SET is_peer_review = ? WHERE id_data = ?", (is_peer_review, id_data))
    elif is_published is not None:
        db.execute("UPDATE projects SET is_published = ?, project_version = ? WHERE id_data = ?", (is_published, id_data))
    if commit:
        db.commit()
    return


def _set_analysis_status(id_analysis, is_peer_review=None, is_published=None, project_version=0, commit=False):
    """Sets the publication bits (is_peer_review, is_published) for the specified analysis. If None, their value is
    unchanged."""
    db = get_db()
    if is_peer_review is not None and is_published is not None:
        db.execute("UPDATE analysis "
                   "SET is_peer_review = ?, is_published = ?, project_version = ? "
                   "WHERE id_analysis = ?", (is_peer_review, is_published, project_version, id_analysis))
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


def _set_hold(id_publication: str, value=1, commit=False):
    """Notes that the project is being held (e.g., waiting for external users)."""
    db = get_db()
    db.execute("UPDATE publications SET hold = ? WHERE id_publication = ?", (value, id_publication))
    if commit:
        db.commit()


def get_publication_author_emails(id_publication: str) -> list:
    """Returns the emails of authors associated with the specified project."""
    db = get_db()
    # Get PSCS users
    pscs_authors = db.execute("SELECT users_auth.email "
                              "FROM users_auth INNER JOIN publications_authors "
                              "ON users_auth.id_user = publications_authors.id_user "
                              "WHERE publications_authors.id_publication = ?", (id_publication,)).fetchall()
    external_authors = db.execute("SELECT email "
                                  "FROM publications_external_authors "
                                  "WHERE id_publication = ?", (id_publication,)).fetchall()
    return [au["email"] for au in pscs_authors + external_authors]


def notify_external_authors(id_publication: str):
    """Notifies external authors that they should input their information."""
    db = get_db()
    external_authors = db.execute("SELECT email "
                                  "FROM publications_external_authors_info "
                                  "WHERE id_publication = ?", (id_publication,)).fetchall()
    addresses = [au["email"] for au in external_authors]
    from pscs.templates.misc import external_author_info_request
    for addr in addresses:
        url_signer = URLSafeTimedSerializer(secret_key=current_app.config["SECRET_KEY"], salt="external_author")
        token = url_signer.dumps(addr)
        author_url = url_for("projects.external_authors_info", id_publication=id_publication, token=token)
        external_author_url = build_full_url(author_url)
        notif_email = external_author_info_request.format(pscs_url=current_app.config["CURRENT_URL"],
                                                          author_info_url=external_author_url)
        subject = "PSCS Project Publication - Authorship Information Required"
        send_email([addr], subject=subject, body=notif_email)
    return


def build_full_url(local_url: str) -> str:
    """Builds the full URL including https://[domain] before the local URL."""
    # If current_url is localhost, use non-secure http
    protocol = ""
    if not current_app.config["CURRENT_URL"].startswith("https://"):
        protocol = "https://"
    if current_app.config["CURRENT_URL"].split(":")[0] == "localhost":
        protocol = "http://"
    return protocol + current_app.config["CURRENT_URL"] + local_url


@bp.route("<id_publication>/external_author/<token>", methods=["GET", "POST"])
def external_authors_info(id_publication, token):
    """
    Page for requesting an external author's information for publishing a project.
    Parameters
    ----------
    id_publication : str
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
    proj = db.execute("SELECT id_publication "
                      "FROM publications_external_authors_info "
                      "WHERE id_publication = ? AND email = ?", (id_publication, external_author_email)).fetchone()
    if proj is None:  # id_project & email don't match user
        flash("The URL is invalid.")
        return redirect(url_for("pscs.index"))
    id_project = db.execute("SELECT id_project FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()["id_project"]
    if request.method == "GET":
        # Get project info to display to user
        project_summary = get_project_summary(db, id_project)
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
        flash_msg = []
        # Form validation
        if not noPHI or not dataUse:
            flash_msg.append("You must confirm that you have read, understood, and agreed to the PHI conditions and Data Use Agreement.")
        if len(author_name) == 0:
            flash_msg.append("You must supply your name.")
        if len(affiliations) == 0:
            flash_msg.append("You must supply an affiliation. If unaffiliated, enter ""Unaffiliated.""")
        if len(flash_msg) != 0:
            flash(" - ".join(flash_msg))
            return redirect(url_for("projects.external_author_info", id_project=id_project, token=token))

        # Collect external author info
        author_info = {"email": escape(external_author_email),
                       "name": author_name,
                       "affiliations": affiliations,
                       "noPHI": noPHI,
                       "dataUse": dataUse,
                       "ip": escape(request.remote_addr)}
        store_external_author_info(id_publication=id_publication, author_information=author_info)
        done = proceed_if_possible(id_project)
        if done:
            flash("Project is now ready.")
            contact_user = db.execute("SELECT id_user_contact FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()["id_user_contact"]
            current_info = _get_publication_info(db, id_project=id_project)[-1]  # only get latest
            peer_review_password = generate_peer_review_password(id_publication=id_publication)
            notify_user_peer_review_password(id_user=contact_user,
                                             publication_info=current_info,
                                             peer_review_password=peer_review_password)
            notify_authors_peer_review(publication_info=current_info,
                                       publishing_id_user=contact_user)
            return redirect(url_for("publications.request_publication", id_publication=id_publication))
        else:
            flash("Project is now waiting for other external authors. An email notification will be sent when the "
                  "project is ready.")
        return redirect(url_for("pscs.index"))
    return


def proceed_if_possible(id_publication: str):
    """Proceeds with publication stage if all external authors have approved."""
    # Check if all external authors have approved
    db = get_db()
    confirmed = db.execute("SELECT confirmed_association "
                           "FROM publications_external_authors_info "
                           "WHERE id_publication = ?", (id_publication,)).fetchall()
    # Checks if any of the authors have not confirmed.
    missing_confirmation = False in list(map(lambda x: x["confirmed_association"], confirmed))
    if missing_confirmation:
        return False  # not everyone has confirmed; exit
    db.execute("UPDATE publications SET hold = FALSE WHERE id_publication = ?", (id_publication,))
    db.commit()
    return True


def store_external_author_info(id_publication: str,
                               author_information: dict):
    """
    Stores the external author information into the DB.
    Parameters
    ----------
    id_publication : str
        ID of the publication for the author.
    author_information : dict
        Dict containing the following fields: email, name, affiliations, noPHI, dataUse

    Returns
    -------
    None
    """
    db = get_db()
    db.execute("UPDATE publications_external_authors_info SET name = ?, confirmed_datause = ?, confirmed_phi = ?, confirmed_association = 1 WHERE id_publication = ? AND email = ?",
               (author_information["name"], author_information["dataUse"], author_information["noPHI"], id_publication, author_information["email"]))
    for idx, aff in enumerate(author_information["affiliations"]):
        db.execute("INSERT INTO publications_external_author_affiliation (id_publication, email, affiliation, affiliation_order) VALUES(?,?,?,?)",
                   (id_publication, author_information["email"], aff, idx))
    db.commit()
    return


def get_project_summary(db: sqlite3.Connection, id_project: str) -> dict:
    """
    Builds a summary of the specified project.
    Parameters
    ----------
    db : sqlite3.Connection
        DB connection
    id_project : str
        ID of the project to summarize.
    Returns
    -------
    dict
        Summary of the project.
    """
    project_summary = dict()
    project_summary |= _get_project_user_summary(db, id_project)
    project_summary |= _get_project_data_summary(db, id_project)
    project_summary |= _get_project_analysis_summary(db, id_project)
    project_summary |= _get_project_jobs_summary(db, id_project)

    user_role = db.execute("SELECT role "
                           "FROM projects_roles "
                           "WHERE id_project = ? AND id_user = ?", (id_project, session["id_user"])).fetchone()
    project_summary["user_role"] = user_role["role"]
    project_summary["publication_status"] = _get_project_publication_status(db, id_project)
    project_summary |= dict(db.execute("SELECT name_project, description "
                                  "FROM projects "
                                  "WHERE id_project = ?", (id_project,)).fetchone())
    return project_summary


def _get_project_user_summary(db, id_project: str) -> dict:
    """Returns a dict containing the entry "users"  of users associated with the specified project."""
    users = db.execute("SELECT users_auth.name_user "
                       "FROM users_auth INNER JOIN projects_roles "
                       "ON users_auth.id_user = projects_roles.id_user "
                       "WHERE projects_roles.id_project = ?", (id_project,)).fetchall()
    summary = {"users": [u["name_user"] for u in users]}
    return summary


def _get_project_data_summary(db, id_project: str) -> dict:
    """Returns the names of the datasets related to the id_project"""
    data = db.execute("SELECT file_name "
                      "FROM data "
                      "WHERE id_project = ?", (id_project,)).fetchall()
    summary = {"data_names": [d["file_name"] for d in data]}
    return summary


def _get_project_data_info(db, id_project: str) -> dict:
    """Returns full info of data related to the specified project."""
    data = db.execute("SELECT data.id_data, users_auth.name_user, data.file_path, data.data_type, data.data_uploaded_time, data.file_hash, data.file_name "
                      "FROM users_auth INNER JOIN data "
                      "ON users_auth.id_user = data.id_user "
                      "WHERE data.id_project = ? "
                      "ORDER BY users_auth.name_user ASC", (id_project,)).fetchall()
    info = {"data": [dict(d) for d in data]}
    for d in info["data"]:
        d["file_path"] = os.path.basename(d["file_path"])
    return info


def _get_project_analysis_summary(db, id_project: str) -> dict:
    """Returns a summary of the analyses related to the id_project"""
    analyses = db.execute("SELECT analysis_name "
                          "FROM analysis "
                          "WHERE id_project = ?", (id_project,)).fetchall()
    summary = {"analysis": [a["analysis_name"] for a in analyses]}
    return summary


def _get_project_analysis_info(db, id_project: str) -> dict:
    """Returns detailed info about analyses associated with the project"""
    analyses = db.execute("SELECT A.analysis_name, A.id_analysis "
                          "FROM analysis AS A "
                          "WHERE id_project = ?", (id_project,)).fetchall()
    return analyses


def _get_project_jobs_summary(db, id_project: str) -> dict:
    """Returns a summary of jobs"""
    jobs = db.execute("SELECT id_job, is_complete "
                      "FROM submitted_jobs "
                      "WHERE id_project = ?", (id_project,)).fetchall()
    summary = {"jobs": [j for j in jobs if not j["is_complete"]]}
    summary["completed_jobs"] = [j for j in jobs if j["is_complete"]]
    return summary


def _get_project_jobs_info(db, id_project: str) -> (dict, dict):
    """Returns the details about computational jobs."""
    running_jobs = db.execute("SELECT J.submitted_resource, J.server_response, J.date_submitted, U.name_user, A.analysis_name "
                              "FROM submitted_jobs AS J INNER JOIN users_auth AS U ON J.id_user = U.id_user "
                              "INNER JOIN analysis AS A ON J.id_analysis = A.id_analysis "
                              "WHERE J.id_project = ? AND J.is_complete = 0 "
                              "ORDER BY J.date_submitted ASC", (id_project,)).fetchall()
    completed_jobs = db.execute("SELECT J.submitted_resource, J.server_response, J.date_submitted, J.id_job, U.name_user, "
                                "A.analysis_name, J.date_completed "
                                "FROM submitted_jobs AS J INNER JOIN users_auth AS U ON J.id_user = U.id_user "
                                "INNER JOIN analysis AS A ON J.id_analysis = A.id_analysis "
                                "WHERE J.id_project = ? AND J.is_complete = 1 "
                                "ORDER BY J.date_submitted ASC", (id_project,)).fetchall()
    completed_jobs_dicts = [dict(c) for c in completed_jobs]
    info = {"jobs": running_jobs}
    info["completed_jobs"] = completed_jobs_dicts
    return info


def _get_project_results_info(db, id_project: str) -> list:
    """Returns details about results associated with the project"""
    results = db.execute("SELECT R.file_path, A.analysis_name, R.is_published, R.is_peer_review, R.result_type, R.is_interactive, R.interactive_tag, R.id_result, R.file_name "
                         "FROM results AS R INNER JOIN analysis AS A ON R.id_analysis = A.id_analysis "
                         "WHERE R.id_project = ?", (id_project,)).fetchall()
    info = [dict(r) for r in results]
    for r in info:
        if r["is_published"]:
            r["publication_status"] = "published"
        elif r["is_peer_review"]:
            r["publication_status"] = "peer review"
        else:
            r["publication_status"] = "private"
    return info


def _get_project_management_info(db, id_project: str) -> dict:
    """Returns info about managing the project"""
    users = db.execute("SELECT PR.role, U.name_user "
                       "FROM projects_roles AS PR INNER JOIN users_auth as U ON PR.id_user = U.id_user "
                       "WHERE PR.id_project = ?", (id_project,)).fetchall()
    papers = db.execute("SELECT doi, url, title, year, author_str FROM project_papers WHERE id_project = ?", (id_project,)).fetchall()
    description = db.execute("SELECT description FROM projects WHERE id_project = ?", (id_project,)).fetchone()["description"]
    project_info = {"users": users, "description": description}
    project_info["publication_status"] = _get_project_publication_status(db, id_project)
    project_info["papers"] = papers
    return project_info


def _get_publication_info(db, id_project: str) -> list:
    publications = db.execute("SELECT P.title, P.description, P.creation_time, P.id_publication, P.status, "
                              " PS.id_publication_short FROM publications AS P INNER JOIN "
                              "publications_shortid AS PS ON P.id_publication = PS.id_publication "
                              "WHERE P.id_project = ? "
                              "ORDER BY P.creation_time ASC", (id_project,)).fetchall()
    pub_dict = [dict(p) for p in publications]
    for d in pub_dict:
        d["url"] = url_for("publications.request_publication", id_publication=d["id_publication"])
        d["shorturl"] = url_for("publications_short.short_publication", id_publication_short=d["id_publication_short"])
    return pub_dict


@bp.route("/<id_project>/logs/<id_job>", methods=["GET"])
@login_required
def get_logs(id_project, id_job):
    if check_user_permission(permission_name="data_read",
                             permission_value=1,
                             id_project=id_project):
        # User has permission to access data
        db = get_db()
        stdout_log, stderr_log = read_logs(db, id_job)
        return jsonify({"stdout": stdout_log, "stderr": stderr_log})


@bp.route("/<id_project>/tabs/<tab>", methods=["GET"])
@login_required
def get_tab_info(id_project, tab):
    if check_user_permission(permission_name="data_read",
                             permission_value=1,
                             id_project=id_project):
        # User has permission to access data
        db = get_db()
        if tab == "data":
            data = _get_project_data_info(db, id_project)
            return render_template("pscs/project_tabs/data.html", data_info=data["data"], id_project=id_project)
        elif tab == "jobs":
            jobs = _get_project_jobs_info(db, id_project)
            return render_template("pscs/project_tabs/jobs.html",
                                   running_jobs=jobs["jobs"],
                                   completed_jobs=jobs["completed_jobs"])
        elif tab == "analysis":
            analysis = _get_project_analysis_info(db, id_project)
            default_analysis = get_default_analysis()
            return render_template("pscs/project_tabs/analysis.html", analysis=analysis, default_analysis=default_analysis)
        elif tab == "results":
            results_info = _get_project_results_info(db, id_project)
            interactive_results_info = dd(list)
            for r in results_info:
                r["file_path"] = r["file_path"][len(current_app.config["INSTANCE_PATH"]):]
                if r["is_interactive"]:
                    interactive_results_info[r["interactive_tag"]].append(r)
            analysis = _get_project_analysis_info(db, id_project)
            return render_template("pscs/project_tabs/results.html", results=results_info, analysis=analysis, interactive_results=interactive_results_info)
        elif tab == "project_management":
            project_info = _get_project_management_info(db, id_project)
            publication_info = _get_publication_info(db, id_project)
            for p in publication_info:
                p["url"] = build_full_url(p["url"])
                p["shorturl"] = build_full_url(p["shorturl"])
            return render_template("pscs/project_tabs/project_management.html",
                                   project_info=project_info,
                                   publication_info=publication_info)
        elif tab == "developer":
            f = open(current_app.config["SSH_PUBKEY"], "r")
            pubkey = f.read()
            f.close()
            return render_template("pscs/project_tabs/developer.html",
                                   ssh_pubkey=pubkey)
    return

@bp.route("/<id_project>/pipeline_import", methods=["POST"])
@login_required
def import_pipeline(id_project):
    """Processes the request for an analysis and imports it to the specified project if allowed."""
    if request.method == "POST":
        # Get pipeline ID from the form
        id_analysis = request.form["id_analysis"]
        if check_analysis_published(id_analysis=id_analysis) or \
            check_user_permission(permission_name="analysis_read",
                                 permission_value=1,
                                 id_project=id_project):
            copy_pipeline(id_analysis, id_project_destination=id_project)
            flash("Pipeline imported.")
        else:
            flash("Either the pipeline doesn't exist or you don't have permission to import it.")
    return redirect(url_for("projects.project", id_project=id_project))

@bp.route("/<id_project>/public_pipeline_import", methods=["POST"])
@login_required
def import_public_pipeline(id_project):
    if request.method == "POST":
        data = request.json
        user_id_project = data["id_project"]
        id_analysis = data["id_analysis"]
        if check_analysis_published(id_analysis=id_analysis) and \
            check_user_permission(permission_name="analysis_write",
                                  permission_value=1,
                                  id_project=user_id_project):
            # Copy!
            copy_pipeline(id_analysis=id_analysis, id_project_destination=user_id_project)
            return jsonify({"success": 1, "message": ""})
    return jsonify({"success": 0, "message": "Insufficient permissions."})


@bp.route("/<id_project>/file_request", methods=["POST"])
@limiter.limit("1/minute")
def file_request(id_project):
    if request.method == "POST":
        # First check that the user is allowed to download these
        req_files = request.json
        db = get_db()
        file_info = []
        for id_data in req_files:
            # Check if data is public
            data = db.execute("SELECT is_published, is_peer_review, file_path, file_name FROM data WHERE id_data = ?", (id_data,)).fetchone()
            data_read = check_user_permission("data_read", permission_value=1, id_project=id_project)
            if data_read or data["is_published"] == 1 or data["is_peer_review"] == 1:
                file_info.append({"file_path": data["file_path"], "file_name": data["file_name"]})
                continue
            else:
                return jsonify({"success": 0, "message": "You do not have permission to access the requested files."}), 403
        else:
            if len(file_info) == 0:
                return jsonify({"success": 0, "message": "The files could not be sent."})
            # Get file paths of each file; zip into archive, send to user
            try:
                # Get project name
                zip_name = zip_files(file_info)
                return send_file(zip_name, as_attachment=True, download_name=f"data_{id_project}.zip",
                                 mimetype="application/zip")
            finally:
                os.remove(zip_name)
    return jsonify({"success": 0, "message": "Unable to send file."})


@bp.route("/<id_project>/results_request", methods=["POST"])
@limiter.limit("1/minute")
def results_request(id_project):
    if request.method == "POST":
        # Get request info
        req_results = request.json
        id_analysis = req_results["id_analysis"]
        db = get_db()
        # Check whether user is allowed to get results
        if is_user_allowed_read_results(id_analysis):
            # Get results for the analysis; zip and send
            # First get id_results
            results_data = db.execute("SELECT id_result, file_path FROM results "
                                    "WHERE id_analysis = ?", (id_analysis,)).fetchall()
            results_list = []
            zip_name = None
            for res in results_data:
                results_list.append({"file_path": res["file_path"]})
            try:
                zip_name = zip_files(results_list)
                return send_file(zip_name, as_attachment=True, download_name=f"results_{id_analysis}.zip",
                                 mimetype="application/zip")
            finally:
                if zip_name is not None:
                    os.remove(zip_name)
    return jsonify({"success": 0, "message": "Permission denied."}), 403


def is_user_allowed_read_results(id_analysis):
    db = get_db()
    id_project = db.execute("SELECT id_project FROM analysis WHERE id_analysis = ?", (id_analysis,)).fetchone()["id_project"]
    pub_status = _get_project_publication_status(db, id_project)
    if check_user_permission(permission_name="data_read",
                             permission_value=1,
                             id_project=id_project) or \
        pub_status == "public" or \
            (pub_status == "peer review" and "project_review" in session.keys() and id_project in session["project_review"]):
        return True
    return False

def zip_files(file_info_list: Collection[dict]):
    """
    Zips several files into a zip archive. Elements should be dicts with "file_path" and "file_name" keys.
    Parameters
    ----------
    file_info_list : Collection[dict]
        List of dicts with "file_path" and "file_name" keys defined. "file_path" points to the location of the file
        on disk. "file_name" is the name that the file should be assigned in the file.

    Returns
    -------
    str
        Path to the zip archive.
    """
    zip_archive = tempfile.NamedTemporaryFile(suffix=".zip", prefix="pscs_filereq_", delete=False)
    zip_writer = zipfile.ZipFile(zip_archive, "w")
    for finfo in file_info_list:
        zip_writer.write(finfo["file_path"], arcname=os.path.basename(finfo["file_path"]))
    zip_writer.close()
    return zip_archive.name


def copy_pipeline(id_analysis, id_project_destination):
    """Imports a pipeline into the specified project."""
    # NOTE: Currently, the node_file and parameter_file are the same file. If that changes, this function
    # needs to be updated.
    db = get_db()
    pipeline = dict(db.execute("SELECT analysis_name, node_file, analysis_hash, is_validated, initial_pscs_version "
                               "FROM analysis WHERE id_analysis = ?", (id_analysis,)).fetchone())
    pipeline_inputs = db.execute("SELECT node_id, node_name FROM analysis_inputs WHERE id_analysis = ?", (id_analysis,)).fetchall()
    # Get new id for analysis
    new_id = get_unique_value_for_field(db, "id_analysis", "analysis")
    node_file = pipeline["node_file"]
    # Copy files to current project
    project_path = current_app.config["PROJECTS_DIRECTORY"].format(id_project=id_project_destination)
    new_node_file = join(project_path, new_id + ".json")
    shutil.copy(node_file, new_node_file)
    if pipeline:
        try:
            db.execute("INSERT INTO analysis (id_analysis, id_project, analysis_name, node_file, parameter_file, "
                       "analysis_hash, is_validated, initial_pscs_version) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                       (new_id, id_project_destination, pipeline["analysis_name"], new_node_file, new_node_file,
                        pipeline["analysis_hash"], pipeline["is_validated"], pipeline["initial_pscs_version"]))
            for p_input in pipeline_inputs:
                # Get new IDs, insert into analysis_inputs
                id_input = get_unique_value_for_field(db, "id_input", "analysis_inputs")
                db.execute("INSERT INTO analysis_inputs (id_input, id_analysis, node_id, node_name) VALUES (?,?,?,?)", (id_input, new_id, p_input["node_id"], p_input["node_name"]))
            db.commit()
        except Exception as e:
            # Clean up file if something went wrong.
            os.remove(new_node_file)
            raise e
    return new_id


def make_shortid(db: sqlite3.Connection,
                 id_publication: str,
                 num_chars: int = 5,
                 charset: Optional[str] = None):
    """
    For a project that is getting published, get a short ID that can be used for sharing links.
    Parameters
    ----------
    db : sqlite3.Connection
        DB against which to check for the short ID.
    id_publication : str
        ID of the publication for which to create the link.
    num_chars : int
        Number of characters for the shortlink. Default: 5.
    charset: Collection[str]
        The character set from which to choose to make the short ID. Default: 0-9, a-z, A-Z.

    Returns
    -------
    bool
        Whether the operation succeeded.
    """
    if charset is None:
        charset = string.digits + string.ascii_letters
    attempt_count = 0
    while attempt_count < 50:
        attempt_count += 1
        short_id = "".join(random.sample(charset, num_chars))
        # Ensure shortid is unique
        count = db.execute("SELECT COUNT(*) FROM publications_shortid WHERE id_publication_short = ?", (short_id,)).fetchone()[0]
        if count == 0:
            db.execute("INSERT INTO publications_shortid (id_publication, id_publication_short) VALUES(?,?)", (id_publication, short_id))
            db.commit()
            return short_id
    raise IOError(f"Unable to create short ID for publication {id_publication}")

@bp.route("/<id_project>/link_papers", methods=["POST"])
def link_papers(id_project):
    """Processes a request to link papers to a project. Papers should be identified by their DOI. An http GET request
    is used to retrieve relevant info."""
    if request.method == "POST":
        # check permission
        if not check_user_permission("project_management", 1, id_project):
            return jsonify({"status": "error", "msg": "You do not have permission to link papers."}), 403
        doi_list = request.json["doi"]
        if len(doi_list) > current_app.config["MAX_PAPERS"]:
            # Too many
            return jsonify({"status": "error", "msg": "Too many DOI values."}), 422
        # validate doi in doi_list
        for doi in doi_list:
            # Validate; if all are valid, then ping doi.org
            if not validate_doi(doi):
                return jsonify({"status": "error", "msg": "Invalid DOI.", "doi": doi}), 422
        # Ping DOI
        url_list = ["https://doi.org/" + doi for doi in doi_list]
        db = get_db()
        # Check how many are currently linked
        num_papers = db.execute("SELECT COUNT(ALL) FROM project_papers WHERE id_project = ?", (id_project,)).fetchone()[0]
        if num_papers + len(doi_list) > current_app.config["MAX_PAPERS"]:
            return jsonify({"status": "error", "msg": "Projects can only have " + current_app.config["MAX_PAPERS"] + " linked papers."}), 422
        for url, doi in zip(url_list, doi_list):
            req = requests.post(url, headers={"Accept": "application/x-bibtex"})
            if req.status_code != 200:
                return jsonify({"status": "error", "msg": "Unable to get DOI: " + doi})
            req_content = req.content.decode("utf-8")
            bibtex_parser = bibtex.Parser()
            data = bibtex_parser.parse_string(req_content)
            for id_bib in data.entries:
                fields = data.entries[id_bib].fields
                persons = data.entries[id_bib].persons["author"]
                num_authors = len(persons)
                author_list = []
                for author in persons[:3]:
                    firsts = "".join([f[0].upper() + ". " for f in author.first_names])
                    middles = "".join([m[0].upper() + ". " for m in author.middle_names])
                    lasts = " ".join(author.last_names)
                    name_str = " ".join([firsts, middles, lasts]) if len(middles) > 0 else " ".join([firsts, lasts])
                    author_list.append(name_str)
                author_str = ", ".join(author_list)
                # there should be only one entry; ignore if there are multiple
                break
            try:
                db.execute("INSERT INTO project_papers (id_project, doi, url, title, year, author_str, num_authors) VALUES (?,?,?,?,?,?,?)",
                           (id_project, escape(doi), url, escape(fields["title"]), escape(fields["year"]), escape(author_str), num_authors))
            except sqlite3.IntegrityError:
                pass
            if len(url_list) > 1:
                time.sleep(0.05)  # avoid spamming the server
        db.commit()
        return jsonify({"status": "success", "msg": ""}), 200


def validate_doi(doi: str):
    """Checks whether the input string is a valid DOI."""
    # Format should be 123123/XYZ456; 10.XXX/12213 is valid, and alphabetical chars can appear in the second part
    if "/" not in doi:
        return False
    if doi.startswith("http"):
        return False
    doi_split = doi.split("/")
    if len(doi_split) != 2:
        return False
    if not doi_split[0].replace(".", "").isnumeric():
        return False

    if "." in doi:
        f = open(join(current_app.config["STATIC_DIRECTORY"], "tlds.txt"), "r")
        tlds = f.read().splitlines()
        for tld in tlds:
            if tld in doi:
                return False
    return True


@bp.route("/<id_project>/remove_paper", methods=["POST"])
def remove_paper(id_project):
    if request.method == "POST":
        if not check_user_permission("project_management", 1, id_project):
            return jsonify({"status": "error", "msg": "You do not have permission to modify this project"}), 403
        doi = request.json["doi"]
        db = get_db()
        db.execute("DELETE FROM project_papers WHERE id_project = ? AND doi = ?", (id_project, doi))
        db.commit()
        return jsonify({"status": "success", "msg": ""}), 200


@bp.route("/<id_project>/update_figure", methods=["POST"])
def update_figure(id_project):
    if request.method == "POST":
        if not check_user_permission("analysis_write", 1, id_project):
            return jsonify({"status": "error", "msg": "You do not have permission to modify this project"}), 403
        id_result = request.json["id_result"]

        # check that there is an associated .fig
        db = get_db()
        file_path_fig = db.execute("SELECT file_path_fig FROM results_figures WHERE id_result = ?", (id_result,)).fetchone()
        # Load & update figure with new parameters
        f = open(file_path_fig["file_path_fig"], "rb")
        fig = pkl.load(f)
        f.close()
        ax = fig.get_axes()[0]
        for k, v in request.json.items():
            if k == "id_result":
                continue
            else:
                print(f"updating {k}")
                update_property(ax, k, v)
        img_path = db.execute("SELECT file_path FROM results WHERE id_result = ?", (id_result,)).fetchone()["file_path"]
        print(f"saving to {img_path}")
        fig.savefig(img_path)
        # Update pkl file, too
        f = open(file_path_fig["file_path_fig"], "wb")
        pkl.dump(fig, f)
        f.close()
        print("SUCCESS")
        return jsonify({"status": "success", "msg": ""}), 200

def update_property(ax, prop, value):
    if value is None:
        return
    if prop == "xlabel":
        ax.set_xlabel(value)
    elif prop == "xlabel_fontsize":
        ax.set_xlabel(ax.get_xlabel(), fontdict={"fontsize": value})
    elif prop == "ylabel":
        ax.set_ylabel(value)
    elif prop == "ylabel_fontsize":
        ax.set_ylabel(ax.get_ylabel(), fontdict={"fontsize": value})
    elif prop == "axhline":
        ax.axhline(y=float(value), color="black")
    elif prop == "axvline":
        ax.axvline(x=float(value), color="black")
    elif prop == "figure_xsize":
        fig = ax.get_figure()
        fig.set_size_inches(float(value), fig.get_size_inches()[1])
    elif prop == "figure_ysize":
        fig = ax.get_figure()
        fig.set_size_inches(fig.get_size_inches()[0], float(value))
    elif prop == "title":
        ax.set_title(value)
    elif prop == "title_fontsize":
        ax.set_title(ax.get_title(), fontdict={"fontsize": value})
    return
