from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, send_from_directory
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename

from scpp.auth import login_required
from scpp.db import get_db
from scpp.metadata.metadata import get_metadata
import os
import uuid
from scpp.analysis import single_cell
from flask import send_from_directory

default_role = 'owner'
bp = Blueprint('scpp', __name__)
rootdir = '/home/lex/flask-tutorial/'
UPLOAD_FOLDER = '/home/lex/flask-tutorial/upload/'
ALLOWED_EXTENSIONS = {'csv', 'tsv'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.join(UPLOAD_FOLDER, "{userid}")
app.config['COMMON_DIRECTORY'] ='/home/lex/flask-tutorial/common/'
app.config['RESULTS_DIRECTORY'] = os.path.join("scpp","static", "results", "{userid}")


@bp.route('/')
def index():
    db = get_db()
    if g.user is not None:
        user_uploads = db.execute(
            'SELECT file_path, sample_count, gene_count, file_hash FROM data WHERE'
            ' userid = ?',
            (g.user['id'],)
        ).fetchall()
        img = 'static/highest_expr_genes_test.jpg'
        return render_template('scpp/index.html', data=user_uploads, img=img)
    else:
        return render_template('scpp/index.html', data=None)


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
            out_path = app.config['UPLOAD_FOLDER'].format(userid=g.user['id'])
            os.makedirs(out_path, exist_ok=True)
            out_file = os.path.join(out_path, filename)
            file.save(out_file)
            if filename.endswith('tsv') or filename.endswith('csv'):
                meta_dict = get_metadata(out_file)
                sample_count, gene_count = meta_dict['table_dimensions']
                hash_value = meta_dict['table_hash']
                db = get_db()
                db.execute(
                    'INSERT INTO data (userid, file_path, sample_count, gene_count, file_hash)'
                    ' VALUES (?,?,?,?,?)',
                    (g.user['id'], os.path.basename(filename), sample_count, gene_count, hash_value))
                db.commit()
            return redirect(url_for('scpp.profile'))
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
    user_projects = db.execute(
        'SELECT'
    )
    user_uploads = db.execute(
        'SELECT file_path, sample_count, gene_count, file_hash FROM data WHERE'
        ' userid = ?',
        (g.user['id'],)).fetchall()
    return render_template('scpp/profile.html', data=user_uploads)

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
                project_row = db.execute("SELECT id FROM projects WHERE id=(?)", (project_id,)).fetchall()
            user_id = g.user['id']
            db.execute("INSERT INTO projects (id, name, description, userid, role) VALUES (?,?,?,?,?)", (project_id, project_name, project_description, user_id, default_role))
            db.commit()
        return redirect(url_for('scpp.index'))

    return render_template("scpp/create.html")

@bp.route('/analysis', methods=['GET', 'POST'])
@login_required
def analysis():
    # Page needs file list, project info
    userid = g.user['id']
    db = get_db()
    # Get files associated with this user
    data = db.execute("SELECT file_path FROM data WHERE userid=(?)", (userid,)).fetchall()
    if request.method == 'GET':
        return render_template(("scpp/analysis.html"), files=data)
    if request.method == 'POST':
        file = os.path.join(app.config['UPLOAD_FOLDER'].format(userid=userid), request.form['analyze'])
        res_dir = app.config['RESULTS_DIRECTORY'].format(userid=userid)
        os.makedirs(res_dir, exist_ok=True)
        single_cell.analyze(file, res_dir, 'test')
        results = os.listdir(res_dir)
        results = [os.path.join(res_dir, r) for r in results]
        results.sort(key=os.path.getmtime)
        results_sans_static = [r.replace('scpp/','') for r in results]

        # Move to new page
        return render_template("scpp/analysis.html", files=data, results=results_sans_static)

@bp.route('/results/<userid>/<filename>', methods=['GET'])
@login_required
def get_file(userid, filename):
    return send_from_directory(app.config['RESULTS_DIRECTORY'].format(userid=userid), filename)

app.add_url_rule('/upload/<name>', endpoint='download_file', build_only=True)