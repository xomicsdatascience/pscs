from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app
)

from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename
from pscs.auth import login_required
from pscs.db import get_db, get_unique_value_for_field
from pscs.metadata.metadata import get_metadata
from pscs.analysis.pipeline import node_parser
import os
import uuid
import numpy as np
from pscs.analysis import single_cell
from pscs.analysis.dispatching import dispatch
from flask import send_from_directory
import plotly.express as px
import plotly
import pandas as pd
import json
import hashlib
import pathlib

default_role = 'owner'
bp = Blueprint('pscs', __name__)
UPLOAD_FOLDER = 'upload/'
ALLOWED_EXTENSIONS = {'csv', 'tsv'}
PATH_KEYWORD = 'path'


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.join(UPLOAD_FOLDER, "{userid}")
app.config['PROJECTS_DIRECTORY'] = os.path.join("pscs", "static", "projects", "{id_project}")
app.config['RESULTS_DIRECTORY'] = os.path.join(app.config['PROJECTS_DIRECTORY'], "results", "{id_analysis}")


@bp.route('/')
def index():
    db = get_db()
    if g.user is not None:
        # Get project meta data to list for user
        user_projects = db.execute('SELECT projects.name_project, projects.description, projects.num_files,'
                                   ' projects.id_project, projects.num_members, projects_roles.role'
                                   ' FROM projects'
                                   ' INNER JOIN  projects_roles'
                                   ' ON projects.id_project=projects_roles.id_project'
                                   ' AND projects.id_user=? AND projects_roles.id_user=?',
                                   (g.user['id_user'], g.user['id_user']))

        return render_template('pscs/index.html', projects=user_projects)
    else:
        return render_template('pscs/index.html', projects=None)


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route(f"/upload/<user>/<name>", methods=['GET'])
@login_required
def download_file(name, user):
    return send_from_directory(app.config['UPLOAD_FOLDER'].format(username=user), name)


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
            filename = secure_filename(file.filename)
            out_path = app.config['UPLOAD_FOLDER'].format(userid=g.user['id_user'])
            os.makedirs(out_path, exist_ok=True)
            out_file = os.path.join(out_path, filename)
            file.save(out_file)
            if filename.endswith('tsv') or filename.endswith('csv'):
                meta_dict = get_metadata(out_file)
                sample_count, gene_count = meta_dict['table_dimensions']
                hash_value = meta_dict['table_hash']
                id_project = session['CURRENT_PROJECT']
                db = get_db()
                data_row = 'temporary'
                while len(data_row) > 0:
                    id_data = str(uuid.uuid4())
                    data_row = db.execute('SELECT file_hash FROM data WHERE id_data = ?', (id_data,)).fetchall()
                db.execute(
                    'INSERT INTO data (id_data, id_user, id_project, file_path, data_type, file_hash)'
                    ' VALUES (?,?,?,?,?,?)',
                    (id_data, g.user['id_user'], id_project, out_file, 'table', hash_value))
                db.execute('UPDATE projects SET num_files = num_files + 1 WHERE id_project = ? AND id_user = ?',
                           (id_project, g.user['id_user']))
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


@bp.route('/profile', methods=['GET'])
@login_required
def profile():
    db = get_db()
    user_uploads = db.execute(
        'SELECT file_path, file_hash FROM data WHERE'
        ' id_user = ?',
        (g.user['id_user'],)).fetchall()
    return render_template('pscs/profile.html', data=user_uploads)


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
            # Insert into projects table
            project_row = 'temp'
            while len(project_row) > 0:
                project_id = str(uuid.uuid4())
                project_row = db.execute("SELECT id_project FROM projects WHERE id_project=(?)", (project_id,)).fetchall()
            db.execute("INSERT INTO projects (id_project, id_user, name_project, description) VALUES (?,?,?,?)",
                       (project_id, g.user['id_user'], project_name,project_description))
            db.execute("INSERT INTO projects_roles (id_project, id_user, role) VALUES (?,?,?)",
                       (project_id, g.user['id_user'], 'admin'))
            db.commit()
            proj_dir = pathlib.Path(app.config['PROJECTS_DIRECTORY'].format(id_project=project_id))
            proj_dir.mkdir(exist_ok=True)
            results_dir = pathlib.Path(os.path.join(app.config['PROJECTS_DIRECTORY'], 'results').format(id_project=project_id))
            results_dir.mkdir(exist_ok=True)

        return redirect(url_for('pscs.index'))
    return render_template("pscs/create.html")


@bp.route('/analysis', methods=['GET', 'POST'])
@login_required
def analysis():
    # Page needs file list, project info
    userid = g.user['id_user']
    if 'CURRENT_PROJECT' in session.keys():
        id_project = session['CURRENT_PROJECT']
    else:
        redirect(url_for('pscs.index'))
    db = get_db()
    # Get files associated with this user
    data = db.execute("SELECT file_path FROM data WHERE id_user=(?) AND id_project = ?", (userid,id_project)).fetchall()
    data_list = []
    for d in data:
        data_list.append(dict(d))
        data_list[-1]['file_path_basename'] = os.path.basename(data_list[-1]['file_path'])
    if request.method == 'GET':
        return render_template(("pscs/analysis.html"), files=data)
    if request.method == 'POST':
        file = os.path.join(app.config['UPLOAD_FOLDER'].format(userid=userid), request.form['analyze'])
        res_dir = app.config['RESULTS_DIRECTORY'].format(userid=userid)
        os.makedirs(res_dir, exist_ok=True)
        ann_data = single_cell.analyze(file, res_dir, 'test')
        results = os.listdir(res_dir)
        results = [os.path.join(res_dir, r) for r in results]
        results.sort(key=os.path.getmtime)

        results = db.execute('SELECT file_path, result_type, title, description FROM results WHERE id_project = ?', (id_project,)).fetchall()
        results_list = []
        for r in results:
            results_list.append(dict(r))
            results_list[-1]['file_load_path'] = r['file_path'].replace('pscs/', '')

        plot_data = ann_data.obs
        plot_data['UMAP1'] = ann_data.obsm['X_umap'][:, 0]
        plot_data['UMAP2'] = ann_data.obsm['X_umap'][:, 1]
        index_name = plot_data.index.name
        plot_data.reset_index(inplace=True)
        fig = px.scatter(plot_data, x='UMAP1',
                         y='UMAP2',
                         size='total_counts',
                         hover_name=index_name,
                         size_max=30,
                         color=ann_data.obs['leiden'],
                         width=800)

        graph_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        # Move to new page
        return render_template("pscs/analysis.html", files=data_list, results=results_list, graph_json=graph_json)
        # return render_template("pscs/analysis.html", files=data_list, results=results, graph_json=graph_json)


@bp.route('/results/<userid>/<filename>', methods=['GET'])
@login_required
def get_file(userid, filename):
    return send_from_directory(app.config['RESULTS_DIRECTORY'].format(userid=userid), filename)


@bp.route('/project/<id_project>', methods=['GET', 'POST'])
@login_required
def project(id_project):
    if request.method == 'GET':
        db = get_db()
        id_user = g.user['id_user']
        role = db.execute('SELECT role FROM projects_roles WHERE id_project = ? and id_user = ?', (id_project, id_user)).fetchall()
        if len(role) > 0:
            # Display files for this project only
            session['CURRENT_PROJECT'] = id_project
            project_name = db.execute('SELECT name_project FROM projects WHERE id_project = ? and id_user = ?', (id_project, id_user)).fetchall()[0]['name_project']
            # Get analyses
            project_analyses = db.execute('SELECT id_analysis, analysis_name FROM analysis WHERE id_project = ?', (id_project,)).fetchall()
            analyses = {}
            analysis_nodes = {}
            for an in project_analyses:
                analyses[an['id_analysis']] = an['analysis_name']
                # Get input nodes with analysis
                project_inputs_db = db.execute('SELECT node_id, node_name FROM analysis_inputs WHERE id_analysis = ?', (an['id_analysis'],)).fetchall()
                project_inps = {}
                for inp in project_inputs_db:
                    project_inps[inp['node_id']] = inp['node_name']
                analysis_nodes[an['id_analysis']] = project_inps
            # Get files associated with project
            project_data = db.execute('SELECT id_data, file_path FROM data WHERE id_project = ? AND id_user = ?', (id_project, id_user)).fetchall()
            files = {}
            for project_file in project_data:
                files[project_file['id_data']] = os.path.basename(project_file['file_path'])
            return render_template("pscs/project.html", project_name=project_name, analyses=analyses, files=files, analysis_nodes=analysis_nodes)
    return redirect(url_for('pscs.index'))


@bp.route('/run_analysis', methods=['POST'])
@login_required
def run_analysis():
    if request.method == 'POST':
        pipeline_specs = request.json
        id_project = session['CURRENT_PROJECT']
        id_analysis = pipeline_specs['id_analysis']
        # Get analysis json
        db = get_db()
        # TODO: need to confirm that current user can access this analysis
        pipeline_json = db.execute('SELECT node_file FROM analysis WHERE id_analysis = ?',
                                   (pipeline_specs['id_analysis'],)).fetchall()[0]['node_file']

        # This next section gets the relevant paths for input and output.
        # Instead, it should create the output paths (as needed), gather the input files, and trigger for them to be
        # sent out to OSG.
        output_dir = app.config['RESULTS_DIRECTORY'].format(id_project=session['CURRENT_PROJECT'],
                                                            id_analysis=pipeline_specs['id_analysis'])
        pathlib.Path(output_dir).mkdir(exist_ok=True)
        file_ids = pipeline_specs['file_paths']
        for node_id, file_id in file_ids.items():
            buf = db.execute('SELECT file_path FROM data WHERE id_data = ?', (file_id,)).fetchall()[0]
            file_ids[node_id] = buf['file_path']
        # Dispatch to OSP
        dispatch(pipeline_json=pipeline_json,
                 file_ids=file_ids,
                 id_project=id_project,
                 id_analysis=id_analysis,
                 resource='osp')
        return redirect(url_for('pscs.project', id_project=session['CURRENT_PROJECT']))


@bp.route('/pipeline', methods=['GET', 'POST'])
def pipeline_designer():
    if request.method == 'GET':
        f = open("pscs/static/node_data.json", 'r')
        j = json.load(f)
        f.close()
        modules = get_unique_values_for_key(j, 'module')
        node_json = convert_dict_to_list(j)
        proj_dests, user_dests = get_save_destinations()
        return render_template("pscs/pipeline.html", node_json=node_json, modules=modules, proj_dests=proj_dests, user_dests=user_dests)
    elif request.method == 'POST':
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
        pipeline_id = str(uuid.uuid4())
        pipeline_row = 'temporary'
        while len(pipeline_row) > 0:
            id_analysis = str(uuid.uuid4())
            pipeline_row = db.execute('SELECT id_analysis FROM analysis WHERE id_analysis = ?', (id_analysis,)).fetchall()

        if is_dest_project:
            pipeline_dir = app.config['PROJECTS_DIRECTORY'].format(id_project=id_project)
        else:
            pipeline_dir = app.config['UPLOAD_FOLDER'].format(userid=g['id_user'])

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


app.add_url_rule('/upload/<name>', endpoint='download_file', build_only=True)


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
