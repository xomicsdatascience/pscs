import functools
from uuid import uuid4
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)

from werkzeug.security import check_password_hash, generate_password_hash
from pscs.db import get_db

bp = Blueprint('auth', __name__, url_prefix='/auth')

@bp.route('/register', methods=('GET', 'POST'))
def register():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        db = get_db()
        error = None

        if not username:
            error = 'Username is required.'
        if not password:
            error = ' '.join(error, 'Password is required.').strip()

        if error is None:
            try:
                # Generate uuid for user
                # there's a SLIM chance that uid is already in the table; this would result in the user receiving
                # an error that is not reproducible and not at all their fault
                # we'll first check that it's not in there
                uid_row = 'starter'
                while len(uid_row) > 0:
                    # This section is basically impossible to test; it's okay if tests don't have coverage here
                    uid = str(uuid4())
                    uid_row = db.execute("SELECT id_user FROM users WHERE id_user=(?)", (uid,)).fetchall()

                db.execute(
                "INSERT INTO users (id_user, name_user, password) VALUES (?, ?, ?)",
                (uid, username, generate_password_hash(password)))
                db.commit()
            except db.IntegrityError:
                error = f"User {username} is already registered."
            else:
                return redirect(url_for("auth.login"))

        flash(error)
    return render_template("auth/register.html")


@bp.route('login', methods=('GET', 'POST'))
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']

        db = get_db()
        error = None
        user = db.execute(
            'SELECT id_user, name_user, password FROM users WHERE name_user = ?', (username,)
        ).fetchone()

        if user is None:
            error = 'Incorrect username.'
        elif not check_password_hash(user['password'], password):
            error = "Incorrect password."

        if error is None:
            session.clear()
            session['id_user'] = user['id_user']
            return redirect(url_for('index'))
        flash(error)
    return render_template('auth/login.html')


@bp.before_app_request
def load_logged_in_user():
    user_id = session.get('id_user')
    if user_id is None:
        g.user = None
    else:
        g.user = get_db().execute(
            "SELECT id_user, name_user FROM users WHERE id_user = ?", (user_id,)
        ).fetchone()


@bp.route('/logout')
def logout():
    session.clear()
    return redirect(url_for('index'))


def login_required(view):
    @functools.wraps(view)
    def wrapped_view(**kwargs):
        if g.user is None:
            return redirect(url_for('auth.login'))
        return view(**kwargs)
    return wrapped_view
