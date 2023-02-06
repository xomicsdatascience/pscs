import functools
from uuid import uuid4
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app
)

from werkzeug.security import check_password_hash, generate_password_hash
from pscs.db import get_db, get_unique_value_for_field
from authtools.validation.registration import validate_username, validate_password, validate_email, validate_recaptcha,\
                                              send_user_confirmation_email
bp = Blueprint('auth', __name__, url_prefix='/auth')

@bp.route('/register', methods=('GET', 'POST'))
def register():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        password_confirm = request.form['passwordConfirm']
        email = request.form['email']
        recaptcha = request.form['g-recaptcha-response']
        user_ip = request.remote_addr
        print(request.form)
        print(user_ip)
        db = get_db()
        try:
            # Generate uuid for user
            uid = get_unique_value_for_field(db, 'id_user', 'users_auth')

            # If recaptcha is not valid, stop bothering the server
            recaptcha_valid, recaptcha_msg = validate_recaptcha(recaptcha, current_app.config['RECAPTCHA_SERVER'], user_ip)
            print(recaptcha_valid)
            if not recaptcha_valid:
                error = recaptcha_msg
                # Error from recaptcha is relatively useless; just tell user that there was a problem
                flash("Error with reCAPTCHA token.")
                return render_template("auth/register.html")

            uname_valid, uname_msg = validate_username(username, db)
            password_valid, password_msg = validate_password(password, password_confirm)
            email_valid, email_msg = validate_email(email, db)
            error = ''
            if not uname_valid:
                error = uname_msg + ';'
            if not password_valid:
                error += password_msg + ';'
            if not email_valid:
                error += email_msg + ';'
            if len(error) > 0:
                error = error[:-1]
                flash(error)
                return render_template("auth/register.html")
            else:
                db.execute(
                "INSERT INTO users_auth (id_user, name_user, password, email, ip) VALUES (?, ?, ?, ?, ?)",
                    (uid, username, generate_password_hash(password), email), user_ip)
                db.commit()

                # Send mail to user to verify.
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
            'SELECT id_user, name_user, password FROM users_auth WHERE name_user = ?', (username,)
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
            "SELECT id_user, name_user FROM users_auth WHERE id_user = ?", (user_id,)
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
