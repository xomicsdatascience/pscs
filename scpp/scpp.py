from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename

from scpp.auth import login_required
from scpp.db import get_db
from scpp.metadata.metadata import get_metadata
import os
bp = Blueprint('scpp', __name__)
UPLOAD_FOLDER = '/home/lex/flask-tutorial/upload/'
ALLOWED_EXTENSIONS = {'csv', 'tsv'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.join(UPLOAD_FOLDER, "{username}")

@bp.route('/')
def index():
    db = get_db()
    if g.user is not None:
        user_uploads = db.execute(
            'SELECT file_path, sample_count, gene_count, file_hash FROM data WHERE'
            ' username = ?',
            (g.user['username'],)
        ).fetchall()
        return render_template('scpp/index.html', data=user_uploads)
    else:
        return render_template('scpp/index.html', data=None)

@bp.route('/create', methods=('GET', 'POST'))
@login_required
def create():
    if request.method == 'POST':
        title = request.form['title']
        body = request.form['body']
        error = None

        if not title:
            error = 'Requires title'

        if error is not None:
            flash(error)
        else:
            db = get_db()
            db.execute(
                'INSERT INTO post (title, body, author_id)'
                ' VALUES (?,?,?)',
                (title, body, g.user['id'])
            )
            db.commit()
            return redirect(url_for('scpp.index'))
    return render_template('scpp/create.html')

def get_post(id, check_author=True):
    post = get_db().execute(
        'SELECT p.id, title, body, created, author_id, username'
        ' FROM post p JOIN user u ON p.author_id = u.id'
        ' WHERE p.id = ?',
        (id,)
    ).fetchone()

    if post is None:
        abort(404, f"Post id {id} does not exist.")

    if check_author and post['author_id'] != g.user['id']:
        abort(403)

    return post

@bp.route('/<int:id>/update', methods=('GET', 'POST'))
@login_required
def update(id):
    post = get_post(id)

    if request.method == 'POST':
        title = request.form['title']
        body = request.form['body']
        error = None

        if not title:
            error = 'Title required.'

        if error is not None:
            flash(error)
        else:
            db = get_db()
            db.execute(
                'UPDATE post SET title = ?, body = ?'
                ' WHERE id = ?',
                (title, body, id)
            )
            db.commit()
            return redirect(url_for('scpp.index'))
    return render_template('scpp/update.html', post=post)


@bp.route('/<int:id>/delete', methods=('POST',))
@login_required
def delete(id):
    get_post(id)
    db = get_db()
    db.execute('DELETE FROM post WHERE id = ?', (id,))
    db.commit()
    return redirect(url_for('scpp.index'))

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

from flask import send_from_directory
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
            out_path = app.config['UPLOAD_FOLDER'].format(username=g.user['username'])
            os.makedirs(out_path, exist_ok=True)
            out_file = os.path.join(out_path, filename)
            file.save(out_file)
            if filename.endswith('tsv'):
                meta_dict = get_metadata(out_file)
                sample_count, gene_count = meta_dict['table_dimensions']
                hash_value = meta_dict['table_hash']
                db = get_db()
                db.execute(
                    'INSERT INTO data (username, file_path, sample_count, gene_count, file_hash)'
                    ' VALUES (?,?,?,?,?)',
                    (g.user['username'], os.path.basename(filename), sample_count, gene_count, hash_value))
                db.commit()
            return redirect(url_for('scpp.index'))
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <input type=file name=file>
      <input type=submit value=Upload>
    </form>
    '''



app.add_url_rule('/upload/<name>', endpoint='download_file', build_only=True)