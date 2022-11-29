from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session
)
# from dash import Dash
# import dash.dcc as dcc
# import dash.html as html
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename

from pscs.auth import login_required
from pscs.db import get_db
from pscs.metadata.metadata import get_metadata
import os
import uuid
import numpy as np
from pscs.analysis import single_cell
from flask import send_from_directory
import plotly.express as px
import plotly
import pandas as pd
import json

default_role = 'owner'
bp = Blueprint('pscs', __name__)
rootdir = '/home/lex/flask-tutorial/'
UPLOAD_FOLDER = '/home/lex/flask-tutorial/upload/'
ALLOWED_EXTENSIONS = {'csv', 'tsv'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.join(UPLOAD_FOLDER, "{userid}")
app.config['COMMON_DIRECTORY'] ='/home/lex/flask-tutorial/common/'
app.config['RESULTS_DIRECTORY'] = os.path.join("pscs","static", "results", "{userid}")


@bp.route('/')
def index():
    db = get_db()
    if g.user is not None:
        # Get project meta data to list for user
        user_projects = db.execute('SELECT name_project, description, role, num_files, id_project, num_members FROM projects WHERE id_user = ?',
                                   (g.user['id_user'],)).fetchall()
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
        print(request.form)
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
            db.execute("INSERT INTO projects (id_project, id_user, name_project, description, role) VALUES (?,?,?,?,?)",
                       (project_id, g.user['id_user'], project_name,project_description,'admin'))
            db.commit()
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
        results_sans_static = [r.replace('pscs/','') for r in results]

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

@bp.route('/project/<id_project>', methods=['GET'])
@login_required
def project(id_project):
    db = get_db()
    id_user = g.user['id_user']
    role = db.execute('SELECT role FROM projects WHERE id_project = ? and id_user = ?', (id_project, id_user)).fetchall()
    if len(role) > 0:
        session['CURRENT_PROJECT'] = id_project
        project_name = db.execute('SELECT name_project FROM projects WHERE id_project = ? and id_user = ?', (id_project, id_user)).fetchall()[0]
        project_rows = db.execute('SELECT file_path, data_uploaded_time FROM data WHERE id_project = ? AND id_user = ?', (id_project, id_user)).fetchall()
        project_dicts = []
        for r in project_rows:
            project_dicts.append(dict(r))
            project_dicts[-1]['file_path_basename'] = os.path.basename(project_dicts[-1]['file_path'])
        return render_template("pscs/project.html", project=project_name, files=project_dicts)
    return redirect(url_for('pscs.index'))

@bp.route('/pipeline', methods=['GET','POST'])
@login_required
def pipeline_designer():
    return render_template("pscs/pipeline.html")


@bp.route('/test', methods=['GET'])
def drag_test():
    return render_template("pscs/tests.html")
app.add_url_rule('/upload/<name>', endpoint='download_file', build_only=True)