from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, session, current_app, jsonify
)
import flask
from werkzeug.wrappers import Response
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.serving import run_simple

from markdown import markdown
from werkzeug.utils import secure_filename
from pscs.auth import login_required, is_logged_in
from pscs.db import (get_db, get_unique_value_for_field, check_user_permission, check_analysis_published,
                     check_data_published)
from pscs.metadata.metadata import get_metadata
import os
from pscs.transfers.dispatching import dispatch, can_user_submit
from pscs.utils.pipeline.validation import validate_pipeline
from flask import send_from_directory
import json
import hashlib
import pathlib
import sqlite3
from warnings import warn
import requests
import time
from flask_socketio import SocketIO
from pscs import socketio

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
    return render_template("pscs/index.html", posts=posts_markdown)


@bp.route("/cellxgene/<id_result>")
@login_required
def visualize(id_result):
    if "cxg" in session:
        session["cxg"]["id_result"] = id_result
    else:
        session["cxg"] = {"id_result": id_result}
    session.modified = True
    start_cxg(id_result)
    return render_template("pscs/cxg.html", id_result=id_result)


@socketio.on("connect")
def on_connect():
    session["cxg"]["sid"] = request.sid
    return


@socketio.on("disconnect")
def on_disconnect():
    r = requests.post(current_app.config["CXG_URL"] + "/cxg/release_port", json={"id_user": session["id_user"]})
    return


def start_cxg(id_result):
    db = get_db()
    id_project = get_project_from_result(db, id_result)
    if not check_user_permission("data_read", 1, id_project):
        flash("You are not authorized to view that page.")
        return redirect("/")
    result_path = get_result_path_from_id(db, id_result)
    resp = requests.post(current_app.config["CXG_URL"] + "/cxg/request_cxg_port", json={"id_user": session["id_user"]},
                         headers={"Content-Type": "application/json"})
    if not resp.json()["success"]:
        flash(resp.json()["message"])
        return redirect("/")
    port = resp.json()["port"]
    # Start CXG container:
    _ = requests.post(current_app.config["CXG_URL"] + "/cxg/cellxgene",
                      json={"id_user": session["id_user"], "filepath": result_path},
                      headers={"Content-Type": "application/json"})
    session["cxg"] = {"port": str(port), "start_time": time.time()}
    session.modified = True
    return

@bp.route("/static/assets/<path:path>", methods=["GET"])
def redirect_static(path):
    print("Redirecting!")
    return proxy(path="static/assets/" + path)

@bp.route('/cxg/', defaults={'path': ''}, methods=["GET", "POST"])
@bp.route('/cxg/<path:path>', methods=["GET", "POST"])
@login_required
def proxy(path):
    if "cxg" not in session:
        flash("Result not correctly identified; please contact PSCS team to report the issue.")
        return redirect("/")
    if "port" not in session["cxg"].keys() or session["cxg"]["port"] is None:
        r = requests.post(current_app["CXG_URL"] + "/cxg/get_cxg_port", json={"id_user": session["id_user"]})
        port = r.json()["port"]
        session["cxg"]["port"] = port
        if port is None:
            flash("Unable to start CellXGene instance.")
            return redirect("/")
    port = session["cxg"]["port"]
    # if request.url == "/cxg/static/"
    url_query = "?".join([path, request.query_string.decode()])
    if url_query == "?":
        url_query = ""
    response = requests.get(f'http://localhost:{port}/{url_query}')
    excluded_headers = ['content-encoding', 'content-length', 'transfer-encoding', 'connection']
    headers = [(name, value) for (name, value) in
               response.raw.headers.items() if name.lower() not in excluded_headers]
    response = Response(response.content, response.status_code, headers)
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add("Access-Control-Allow-Credentials", "true")
    return response


def get_project_from_result(db, id_result):
    project = db.execute("SELECT id_project FROM results WHERE id_result = ?", (id_result,)).fetchone()
    return project["id_project"]


def get_result_path_from_id(db, id_result):
    result = db.execute('SELECT file_path FROM results WHERE id_result = ?', (id_result,)).fetchone()
    return result["file_path"]


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route(f"/upload/<id_project>/<name>", methods=['GET'])
@login_required
def download_file(id_project, name):
    return send_from_directory(current_app.config['DATA_DIRECTORY'].format(id_project=id_project), name)


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
                return redirect(url_for("projects.project", id_project=session["CURRENT_PROJECT"]))
            filename = secure_filename(file.filename)
            out_path = current_app.config['DATA_DIRECTORY'].format(id_project=session["CURRENT_PROJECT"])
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
            return redirect(url_for('projects.project', id_project=session['CURRENT_PROJECT']))
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <p>Data should be a .csv with each row as an observation (sample), and each column as a feature (gene, protein).</p>
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
                                   "WHERE projects_roles.id_user=?", (g.user['id_user'],))
        invitations_sent = get_project_invitations_sent(g.user["id_user"])
        invitations_received = get_project_invitations_received(g.user["id_user"])
        return render_template('pscs/projects_summary.html', projects=user_projects,
                               invitations_sent=invitations_sent,
                               invitations_received=invitations_received)


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
            data_dir = pathlib.Path(current_app.config["DATA_DIRECTORY"].format(id_project=id_project))
            data_dir.mkdir(exist_ok=True)
            results_dir = pathlib.Path(os.path.join(current_app.config['PROJECTS_DIRECTORY'], 'results').format(id_project=id_project))
            results_dir.mkdir(exist_ok=True)
            return redirect(url_for('projects.project', id_project=id_project))
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
            return_url = url_for('projects.project', id_project=session['CURRENT_PROJECT'])
            return jsonify({"url": return_url, "submit_status": reason, "submit_success": 0})

        if check_analysis_published(id_analysis) or check_user_permission("analysis_execute", 1, id_project):
            pipeline_json = db.execute('SELECT node_file FROM analysis WHERE id_analysis = ?',
                                       (pipeline_specs['id_analysis'],)).fetchall()[0]['node_file']
        else:
            flash("You do not have permission to run the requested analysis.")
            return redirect(url_for("projects.project", id_project=id_project))

        # This next section gets the relevant paths for input and output.
        # Instead, it should create the output paths (as needed), gather the input files, and trigger for them to be
        # sent out to OSG.

        file_ids = pipeline_specs['file_paths']
        file_info = dict()
        for node_id, file_id in file_ids.items():
            if check_data_published(file_id) or check_user_permission("data_read", 1, id_project):
                buf = db.execute('SELECT file_path FROM data WHERE id_data = ?', (file_id,)).fetchall()[0]
                file_info[node_id] = dict()
                file_info[node_id]["file_path"] = buf["file_path"]
                file_info[node_id]["id"] = file_id
            else:
                flash("You do not have permission to use the some of the specified data.")
                return redirect(url_for("projects.project", id_project=id_project))

        # Dispatch to OSP
        pscs_job_id = dispatch(pipeline_json=pipeline_json,
                               id_user=g.user['id_user'],
                               file_info=file_info,
                               id_project=id_project,
                               id_analysis=id_analysis,
                               resource='osp')
        output_dir = current_app.config['RESULTS_DIRECTORY'].format(id_project=session['CURRENT_PROJECT'],
                                                                    id_analysis=pipeline_specs['id_analysis'],
                                                                    id_job=pscs_job_id)
        pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)
        log_dir = current_app.config["LOG_DIRECTORY"].format(id_project=session["CURRENT_PROJECT"],
                                                             id_analysis=pipeline_specs["id_analysis"],
                                                             id_job=pscs_job_id)
        pathlib.Path(log_dir).mkdir(exist_ok=True, parents=True)

        response_json = {"url": url_for('projects.project', id_project=session['CURRENT_PROJECT']),
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


@bp.route("/pipeline/load_analysis", methods=["POST"])
def load_analysis():
    if "loadAnalysis" in request.json:
        # get id of related project to check user permissions
        db = get_db()
        id_analysis = request.json["loadAnalysis"]
        id_project = db.execute("SELECT id_project "
                                "FROM analysis WHERE id_analysis = ?", (id_analysis,)).fetchone()["id_project"]
        has_perm = check_user_permission(permission_name="analysis_read",
                                         permission_value=1,
                                         id_project=id_project)
        if has_perm:
            # has permission to read; go get analysis file and return JSON
            return load_analysis_from_id(id_analysis)
        return {"": ""}


def load_analysis_from_id(id_analysis):
    db = get_db()
    node_file = db.execute("SELECT node_file FROM analysis WHERE id_analysis = ?", (id_analysis,)).fetchone()[
        'node_file']
    f = open(node_file, 'r')
    node_data = json.load(f)
    f.close()
    return node_data


@bp.route('/pipeline', methods=['GET', 'POST'])
@login_required
def pipeline_designer():
    if request.method == 'GET':
        proj_dests, user_dests = get_save_destinations()

        # Fetch analysis that user is allowed to read
        db = get_db()
        analyses = db.execute("SELECT analysis.id_analysis, analysis.analysis_name "
                              "FROM analysis INNER JOIN projects_roles "
                              "ON analysis.id_project=projects_roles.id_project "
                              "WHERE projects_roles.id_user = ? AND projects_roles.analysis_read = 1",
                              (g.user['id_user'],)).fetchall()
        if "CURRENT_PROJECT" not in session.keys():
            return render_template("pscs/pipeline.html", proj_dests=proj_dests,
                                   user_dests=user_dests, analyses=analyses)
        else:
            return render_template("pscs/pipeline.html", proj_dests=proj_dests,user_dests=user_dests, analyses=analyses,
                                   current_project=url_for("projects.project", id_project=session["CURRENT_PROJECT"]))
    elif request.method == 'POST':
        pipeline_summary = request.json
        id_project = pipeline_summary['saveId']
        # First check that the user is allowed to save analyses
        if not check_user_permission("analysis_write", 1, id_project):
            flash("You do not have permission to save analyses to this project.")
            return jsonify({"error": "User does not have permission to save pipeline to this project"}), 403
        # Validate the pipeline
        is_valid, invalid_pipeline_reasons, invalid_node_reasons = validate_pipeline(pipeline_summary["nodes"])
        if not is_valid:
            flash("The pipeline has problems that could prevent it from executing correctly.")

        # Identify input nodes
        input_nodes = {}
        output_nodes = []
        for n in pipeline_summary["nodes"]:
            if n["pscsType"] == "input":
                input_nodes[n["nodeId"]] = n["labelText"]  # labelText is what is displayed to the user

        # Create new pipeline ID
        db = get_db()
        id_analysis = get_unique_value_for_field(db, "id_analysis", "analysis")
        pipeline_id = id_analysis

        pipeline_dir = current_app.config['PROJECTS_DIRECTORY'].format(id_project=id_project)

        output_name = pipeline_id + '.json'
        pipeline_file = os.path.join(pipeline_dir, output_name)
        f = open(pipeline_file, 'w')
        json.dump(pipeline_summary, f, indent=2)
        f.close()

        pipeline_hash = calc_hash_of_file(pipeline_file)
        pipe_name = secure_filename(pipeline_summary['name'])
        # send to database
        db.execute(
            'INSERT INTO analysis '
            '(id_analysis, id_project, analysis_name, node_file, parameter_file, analysis_hash, is_validated)'
            ' VALUES (?,?,?,?,?,?,?)',
            (id_analysis, id_project, pipe_name, pipeline_file, pipeline_file, pipeline_hash, is_valid))
        for inp_id, inp_name in input_nodes.items():
            # get uuid
            id_input = get_unique_value_for_field(db, 'id_input', 'analysis_inputs')
            db.execute(
                'INSERT INTO analysis_inputs (id_input, id_analysis, node_id, node_name)'
                ' VALUES (?,?,?,?)', (id_input, id_analysis, inp_id, inp_name)
            )
        db.commit()
        return render_template("pscs/pipeline.html")


@bp.route("/pipeline/fetch_nodes", methods=["POST"])
def fetch_nodes():
    # This function loads the node data and returns it as a JSON object
    if request.method == "POST":
        package_json = {"packages": []}
        node_path = current_app.config["NODE_DIRECTORY"]
        for node_file in os.listdir(node_path):
            f = open(os.path.join(node_path, node_file), 'r')
            package_json["packages"].append(json.load(f))
            f.close()
        return jsonify(package_json)


def load_node_specs(json_file: str) -> dict:
    """
    Loads the node specs from the specified file. Format is expected to be:
     {pkg_name: {node_name: {parameters: {param_name: [type, default value] }}}}
    Parameters
    ----------
    json_file : str
        Path to the JSON file containing the node specifications.

    Returns
    -------
    dict
        Nested dict with node information.

    Raises
    ------
    ValueError
        If json_file has more than one top-level key.
    """
    f = open(json_file, "r")
    node_data = json.load(f)
    f.close()
    if len(node_data.keys()) > 1:
        raise ValueError("Node file incorrectly-formatted; top-level key should be name of package.")
    return node_data


def load_many_node_specs(json_files: list) -> dict:
    """
    Loads several node specs from a list.
    Parameters
    ----------
    json_files : list
        List of JSON files containing node specs for different packages.

    Returns
    -------
    dict
        Nested dict with node information.

    Warnings
    --------
    UserWarning
         If two or more of the node files have the same package name, the names will be mangled and the function will
         warn the user.
    """
    node_data = {}
    for node_file in json_files:
        node = load_node_specs(node_file)
        pkg_name = list(node.keys())[0]
        old_pkg_name = pkg_name  # in case of name mangling
        if pkg_name in node_data.keys():
            idx = 1
            old_pkg_name = pkg_name
            while f"{old_pkg_name}_{idx}" in node_data.keys():
                idx += 1
            pkg_name = f"{old_pkg_name}_{idx}"
            warn(UserWarning(f"More than one node package have the same name; renaming {old_pkg_name} to {pkg_name}"))
        node_data[pkg_name] = node[old_pkg_name]
    return node_data


@bp.route("/about", methods=["GET"])
def about():
    if request.method == "GET":
        return render_template("about/about.html")


def get_save_destinations():
    # Get all projects that this user belongs to
    db = get_db()
    projs = db.execute("SELECT projects.name_project, projects.id_project "
                       "FROM projects INNER JOIN projects_roles "
                       "ON projects.id_project=projects_roles.id_project "
                       "WHERE projects_roles.id_user=? "
                       "AND projects_roles.analysis_write=1", (g.user["id_user"],)).fetchall()
    proj_dests = {}
    for p in projs:
        proj_dests[p['name_project']] = p['id_project']
    user_dests = {'user': g.user['name_user'], 'id': g.user['id_user']}
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


@bp.route('/projects/<id_project>/results/<id_analysis>/<id_job>/<path:filename>', methods=['GET'])
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
        return public_results(filename, id_project, id_analysis, id_job=id_job)
    elif public_status["status"] == "peer review":
        return review_results(filename, id_project, id_analysis)


# @login_required
def private_results(filename, id_project, id_analysis, id_job):
    if is_logged_in():
        has_perm = check_user_permission("data_read", 1, id_project)
        if not has_perm:
            return
        res_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project,
                                                                 id_analysis=id_analysis,
                                                                 id_job=id_job)
        return send_from_directory(res_dir, secure_filename(filename))
    else:
        return


def review_results(filename, id_project, id_analysis):
    if "project_review" in session.keys() and id_project in session["project_review"]:
        res_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis)
        return send_from_directory(res_dir, secure_filename(filename))
    return


def public_results(filename, id_project, id_analysis, id_job):
    res_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis, id_job=id_job)
    return send_from_directory(res_dir, secure_filename(filename))


def get_project_invitations_received(id_user: str):
    """Returns a list of invitations where the specified user has been invited."""
    db = get_db()
    ids = db.execute("SELECT id_invitation, id_inviter, id_project, time_sent "
                     "FROM project_invitations "
                     "WHERE id_invitee = ?", (id_user,)).fetchall()
    # ids not very useful; convert to username + project title
    invite = []
    for i in ids:
        # get inviter username
        inv = dict(i).copy()
        inv["inviter_name"] = db.execute("SELECT name_user FROM users_auth WHERE id_user = ?", (i["id_inviter"],)).fetchone()["name_user"]
        # get project name
        proj_name = db.execute("SELECT name_project FROM projects WHERE id_project = ?", (i["id_project"],)).fetchone()["name_project"]
        if len(proj_name) > 20:  # in case project name is too long
            proj_name = proj_name[:20] + "[...]"
        inv["name_project"] = proj_name
        invite.append(inv)
    return invite


def get_project_invitations_sent(id_user: str):
    """Returns a list of invitations where the specified user sent the invitation."""
    db = get_db()
    ids = db.execute("SELECT id_invitation, id_invitee, id_project, time_sent "
                      "FROM project_invitations "
                      "WHERE id_inviter = ?", (id_user,)).fetchall()
    invite = []
    for i in ids:
        inv = dict(i).copy()
        inv["invitee_name"] = db.execute("SELECT name_user FROM users_auth WHERE id_user = ?", (i["id_invitee"],)).fetchone()["name_user"]
        proj_name = db.execute("SELECT name_project FROM projects WHERE id_project = ?", (i["id_project"],)).fetchone()["name_project"]
        inv["name_project"] = proj_name
        invite.append(inv)
    return invite

