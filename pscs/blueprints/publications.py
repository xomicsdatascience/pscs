import sqlite3

from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app, send_from_directory,
    jsonify, send_file
)
from pscs.auth import login_required, is_logged_in
from pscs.db import get_db
from werkzeug.security import check_password_hash
from pscs.extensions.limiter import limiter
import os
import zipfile
import tempfile
from typing import Collection
from pscs.pscs import load_analysis_from_id, check_user_permission
from collections import defaultdict as dd
bp = Blueprint("publications", __name__, url_prefix="/publications")
bp_short = Blueprint("publications_short", __name__, url_prefix="/p/")

@bp_short.route("/<id_publication_short>", methods=["GET"])
def short_publication(id_publication_short):
    db = get_db()
    id_publication = db.execute("SELECT id_publication FROM publications_shortid WHERE id_publication_short = ?", (id_publication_short,)).fetchone()["id_publication"]
    return redirect(url_for("publications.request_publication", id_publication=id_publication))



@bp.route("/<id_publication>", methods=['GET', 'POST'])
def request_publication(id_publication):
    if request.method == "GET":
        db = get_db()
        # Check if public or peer review
        status = db.execute("SELECT status FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()
        if status is not None:
            has_peer_reviewed = (status["status"] == "peer review") and ("publication_review" in session.keys() and id_publication in session["publication_review"])
            if status["status"] == "public" or has_peer_reviewed:
                return display_publication(id_publication)  # publication is public; show it
            elif status["status"] == "peer review":
                return render_template("auth/prompt.html", prompt_label="Peer Review Password")
        else:
            flash("There was a problem fetching the requested publication", "error")
            return redirect(url_for("pscs.index"))

    elif request.method == "POST":
        # This only gets executed if the publication is under peer review.
        submitted_password = request.form["password"]
        # Get passhash from db
        db = get_db()
        passhash = db.execute("SELECT peer_password "
                              "FROM publications_peer_review WHERE id_publication = ?", (id_publication,)).fetchone()["peer_password"]
        if check_password_hash(passhash, submitted_password):
            # Password is correct; allow user to access publication without password
            set_session_key("publication_review", id_publication)
            return display_publication(id_publication)
        else:
            flash("Password is incorrect", "error")
            return redirect(url_for("publications.request_publication", id_publication=id_publication))


@bp.route("/<id_publication>/file_request", methods=["POST"])
@limiter.limit("3/minute")
def file_request(id_publication):
    if request.method == "POST":
        # First check that the user is allowed to download these
        req_files = request.json
        db = get_db()
        file_info = []
        for id_data in req_files:
            # Check if data is public
            data_info = get_data_info(db, id_publication, id_data)
            is_public = data_info["status"] == "public"
            is_peer = data_info["status"] == "peer review" and is_peer_allowed(db, id_publication)
            if is_public or is_peer:
                file_info.append({"file_path": data_info["file_path"], "file_name": data_info["file_name"]})
            else:
                return jsonify({"success": 0, "message": "You do not have permission to access the requested files."}), 403
        else:
            # Get file paths of each file; zip into archive, send to user
            try:
                # Get project name
                zip_name = zip_files(file_info)
                return send_file(zip_name, as_attachment=True, download_name=f"data_{id_publication}.zip",
                                 mimetype="application/zip")
            finally:
                os.remove(zip_name)
    return jsonify({"success": 0, "message": "Unable to send file."})


@bp.route("/<id_publication>/file_request_original", methods=["POST"])
@limiter.limit("3/minute")
def file_request_original(id_publication):
    if request.method == "POST":
        # First check that the user is allowed to download these
        req_files = request.json
        db = get_db()
        file_info = []
        for id_data in req_files:
            # Check if data is public
            data_info = get_data_original_info(db, id_publication, id_data)
            is_public = data_info["status"] == "public"
            is_peer = data_info["status"] == "peer review" and is_peer_allowed(db, id_publication)
            if is_public or is_peer:
                file_info.append({"file_path": data_info["file_path"], "file_name": data_info["file_name"]})
            else:
                return jsonify({"success": 0, "message": "You do not have permission to access the requested files."}), 403
        else:
            # Get file paths of each file; zip into archive, send to user
            try:
                # Get project name
                zip_name = zip_files(file_info)
                return send_file(zip_name, as_attachment=True, download_name=f"data_{id_publication}.zip",
                                 mimetype="application/zip")
            finally:
                os.remove(zip_name)
    return jsonify({"success": 0, "message": "Unable to send file."})


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
        zip_writer.write(finfo["file_path"], arcname=finfo["file_name"])
    zip_writer.close()
    return zip_archive.name


def is_peer_allowed(db, id_publication):
    """Checks whether the specified publication can be viewed by the user."""
    # Check if id_publication is in peer review and that the user has previously provided the password
    status = db.execute("SELECT status FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()["status"]
    if status == "peer review" and "publication_review" in session.keys() and id_publication in session["publication_review"]:
        return True
    return False


def get_data_info(db, id_publication, id_data):
    return db.execute("SELECT P.status, D.file_path, D.file_name FROM publications AS P "
                      "INNER JOIN publications_data AS PD ON P.id_publication = PD.id_publication "
                      "INNER JOIN data AS D ON PD.id_data = D.id_data "
                      "WHERE PD.id_data = ? AND P.id_publication = ?", (id_data, id_publication)).fetchone()

def get_data_original_info(db, id_publication, id_data):
    return db.execute("SELECT P.status, D.file_path, D.file_name FROM publications AS P "
                      "INNER JOIN publications_data_original AS PD ON P.id_publication = PD.id_publication "
                      "INNER JOIN data_original AS D ON PD.id_data_original = D.id_data "
                      "WHERE PD.id_data_original = ? AND P.id_publication = ?", (id_data, id_publication)).fetchone()


def display_publication(id_publication):
    """Prepares and returns the requested publication. This function assumes that the user is allowed to access the
    specified publication."""
    db = get_db()
    # Need to get for pub: title, description, authors, author affiliation, data, analyses, results
    # Need to get for user: user projects
    publication_summary = {}
    descriptor = db.execute("SELECT title, description, id_project FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()

    publication_summary |= dict(descriptor)
    publication_summary["authors"] = _get_publication_author_info(db, id_publication)
    publication_summary["affiliations"] = {}
    for author in publication_summary["authors"]:
        publication_summary["affiliations"][author["author_position"]] = author["affiliations"]

    publication_summary["data"] = _get_publication_data(db, id_publication)
    publication_summary["data_original"] = _get_publication_data_original(db, id_publication)
    publication_summary["analysis"] = _get_publication_analysis(db, id_publication)
    results = _get_publication_results(db, id_publication)
    publication_summary["results"] = dd(list)

    for r in results:
        publication_summary["results"][r["id_analysis"]].append(r)

    id_project = db.execute("SELECT id_project FROM publications WHERE id_publication = ?", (id_publication,)).fetchone()["id_project"]
    other_versions = db.execute("SELECT title, description, creation_time, id_publication FROM publications WHERE id_project = ? "
                                "ORDER BY creation_time DESC", (id_project,)).fetchall()
    other_versions = [dict(v) for v in other_versions]
    for v in other_versions:
        v["publication_url"] = url_for("publications.request_publication", id_publication=v["id_publication"])

    publication_summary["versions"] = other_versions
    publication_summary["papers"] = db.execute("SELECT doi, url, title, year, author_str FROM project_papers WHERE id_project = ?",
                        (id_project,)).fetchall()

    # Get user projects for imports
    user_projects = None
    if "id_user" in session and session["id_user"] is not None:
        user_projects = db.execute("SELECT P.id_project, P.name_project "
                                   "FROM projects AS P INNER JOIN "
                                   "projects_roles AS PR ON P.id_project = PR.id_project WHERE "
                                   "PR.id_user = ? AND PR.analysis_write = 1", (session["id_user"],)).fetchall()

    return render_template("pscs/project_public.html",
                           project_summary=publication_summary,
                           user_projects=user_projects)


@bp.route('/<id_publication>/load_analysis', methods=["POST"])
def load_analysis(id_publication):
    if request.method == "POST":
        if "loadAnalysis" in request.json:
            id_analysis = request.json["loadAnalysis"]
            # Verify that the project is indeed public
            db = get_db()
            status = db.execute("SELECT P.status FROM publications AS P "
                                "INNER JOIN publications_analysis AS PA ON P.id_publication = PA.id_publication "
                                "WHERE PA.id_analysis = ? AND P.id_publication = ?", (id_analysis, id_publication)).fetchone()["status"]
            has_peer_review_password = status == "peer review" and "publication_review" in session.keys() and id_publication in session["publication_review"]
            if status == "public" or has_peer_review_password:
                return load_analysis_from_id(id_analysis)
            return {"": ""}

@bp.route("/<id_publication>/request_results", methods=["POST"])
@limiter.limit("1/minute")
def request_analysis_results(id_publication):
    if request.method == "POST":
        id_analysis = request.json["id_analysis"]
        # Check if the analysis is published
        db = get_db()
        publication_status = db.execute("SELECT status FROM publications AS P INNER JOIN publications_analysis AS PA "
                                        "ON P.id_publication = PA.id_publication WHERE PA.id_publication = ?", (id_publication,)).fetchone()["status"]
        can_peer_review = publication_status == "peer review" and "publication_review" in session.keys() and id_publication in session["publication_review"]
        if publication_status == "public" or can_peer_review:
            return send_results_for_analysis(db, id_analysis)
        flash("You do not have permission to download the requested results.", "warning")
        return jsonify({"success": 0, "message": "Permission denied."}), 403


def send_results_for_analysis(db, id_analysis):
    results_data = db.execute("SELECT file_path FROM results WHERE id_analysis = ?", (id_analysis,)).fetchall()
    zip_name = None
    results_list = [{"file_path": r["file_path"], "file_name": os.path.basename(r["file_path"])} for r in results_data]
    try:
        zip_name = zip_files(results_list)
        return send_file(zip_name, as_attachment=True, download_name=f"results_{id_analysis}.zip",
                         mimetype="application/zip")
    finally:
        if zip_name is not None:
            os.remove(zip_name)
    return jsonify({"success": 0, "message": "Permission denied."}), 403


def _get_publication_author_info(db, id_publication) -> list[dict]:
    """Fetches the username, name, author position, and affiliation of all authors associated with the publication"""
    # Get: username, name, affiliation
    # First get users
    authors = db.execute("SELECT USER.id_user, USER.name_user, USER.name, AU.author_position "
                         "FROM users_auth AS USER INNER JOIN publications_authors AS AU "
                         "ON USER.id_user = AU.id_user WHERE AU.id_publication = ? ORDER BY AU.author_position ASC", (id_publication,)).fetchall()
    if authors is None:
        return []
    authors = [dict(a) for a in authors]
    for a in authors:
        # Get affiliations and make string presentable
        aff_list = db.execute("SELECT AFF.affiliation FROM publications_authors_affiliation as AFF "
                                       "WHERE id_user = ? AND id_publication = ? ORDER BY AFF.affiliation_order ASC", (a["id_user"], id_publication)).fetchall()
        a["affiliations"] = ", ".join([aff["affiliation"] for aff in aff_list])
    # Repeat for external authors
    external_authors = db.execute("SELECT email, name, author_position FROM publications_external_authors_info AS EXT "
                                  "WHERE id_publication = ?", (id_publication,)).fetchall()
    ext_authors = [dict(e) for e in external_authors]
    for e in ext_authors:
        aff_list = db.execute("SELECT affiliation FROM publications_external_author_affiliation WHERE id_publication = ? and email = ? ORDER BY affiliation_order ASC", (id_publication, e["email"])).fetchall()
        e["affiliations"] = ", ".join([aff["affiliation"] for aff in aff_list])
    all_authors = authors + ext_authors
    # Reorder list for correct authorship position
    all_authors.sort(key=lambda x: x["author_position"])
    return all_authors


def _get_publication_data(db, id_publication) -> list[dict]:
    """Gets the id, filepath, data type, and file name for the data associated with the publication."""
    data = db.execute("SELECT D.id_data, D.file_path, D.data_type, D.file_name FROM publications_data as PD "
                      "INNER JOIN data AS D ON PD.id_data = D.id_data "
                      "WHERE PD.id_publication = ?", (id_publication,)).fetchall()
    if data is None:
        return []
    return [dict(d) for d in data]

def _get_publication_data_original(db, id_publication) -> list[dict]:
    """Gets the id, filepath, data type, and file name for the original data associated with the publication."""
    data = db.execute("SELECT D.id_data, D.file_path, D.file_name "
                      "FROM data_original as D INNER JOIN publications_data_original as PD ON D.id_data = PD.id_data_original "
                      "WHERE PD.id_publication = ?", (id_publication,)).fetchall()
    print([dict(d) for d in data])
    if data is None:
        return []
    return [dict(d) for d in data]

def _get_publication_analysis(db, id_publication) -> list[dict]:
    """Gets the name and ID of analysis"""
    analysis = db.execute("SELECT A.id_analysis, A.analysis_name FROM publications_analysis as PA "
                          "INNER JOIN analysis AS A ON A.id_analysis = PA.id_analysis"
                          " WHERE PA.id_publication = ?", (id_publication,)).fetchall()
    if analysis is None:
        return []
    return [dict(a) for a in analysis]


def _get_publication_results(db, id_publication) -> list[dict]:
    """Gets the id, analysis, result type, title, description, and interactivity of results associated with
    the publication."""
    result = db.execute("SELECT R.id_result, R.id_analysis, R.description, R.title, R.is_interactive, R.file_path, R.file_name, R.result_type FROM "
                        "publications_results AS PR INNER JOIN results AS R ON R.id_result = PR.id_result "
                        "WHERE PR.id_publication = ?", (id_publication,)).fetchall()
    if result is None:
        return []
    dresult = [dict(r) for r in result]
    for d in dresult:
        if d["file_path"].startswith(current_app.config["INSTANCE_PATH"]):
            d["file_path"] = d["file_path"][len(current_app.config["INSTANCE_PATH"]):]
        if not d["file_path"].startswith(os.sep):
            d["file_path"] = os.sep + d["file_path"]
    return dresult

def set_session_key(session_key: str,
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
    session[session_key] = [value]
    return