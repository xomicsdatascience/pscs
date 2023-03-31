import functools
from uuid import uuid4
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app
)

from werkzeug.security import check_password_hash, generate_password_hash
from pscs.db import get_db, get_unique_value_for_field
from pscs.authtools.authtools.validation.registration import validate_username, validate_password, validate_email, \
                                                             validate_recaptcha, send_user_confirmation_email, decode_token

bp = Blueprint('auth', __name__, url_prefix='/auth')


@bp.route('/register', methods=('GET', 'POST'))
def register():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        password_confirm = request.form['passwordConfirm']
        email = request.form['email']

        user_ip = request.remote_addr
        db = get_db()
        try:
            # Generate uuid for user
            id_user = get_unique_value_for_field(db, 'id_user', 'users_auth')

            if current_app.config["RECAPTCHA_ENABLED"]:  # check that recaptcha is enabled
                recaptcha = request.form['g-recaptcha-response']
                recaptcha_valid, recaptcha_msg = validate_recaptcha(recaptcha, current_app.config['RECAPTCHA_SERVER'], user_ip)
                # If recaptcha is not valid, stop bothering the server
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
                    (id_user, username, generate_password_hash(password), email, user_ip))
                db.commit()

                # Send mail to user to verify.
                if current_app.config["REGISTRATION_REQUIRES_CONFIRMATION"]:
                    send_user_confirmation_email(id_user, email, username)
                    flash(f"Email confirmation sent to {email}; please click the verification link.")
                else:  # email confirmation disabled
                    confirm_user(id_user)
                user_login(username, password)
                # return render_template("pscs/index.html")
                return redirect(url_for("pscs.index"))

        except db.IntegrityError:
            error = f"User {username} is already registered."
        else:
            return redirect(url_for("index"))

        flash(error)
    return render_template("auth/register.html", recaptcha_enabled=current_app.config["RECAPTCHA_ENABLED"])


@bp.route('/confirmation/<token>')
def user_confirmation(token):
    valid_token, token_data, invalid_reason = decode_token(token,
                                                           "confirmation",
                                                           current_app.config["CONFIRMATION_TIMEOUT_SECONDS"])
    if not valid_token:  # Invalid token; exit
        if len(invalid_reason) > 0:
            flash(invalid_reason)
    else:
        # token_data contains id_user; check whether user exists & needs confirming
        db = get_db()
        user_needs_confirmation = db.execute("SELECT confirmed, name_user "
                                             "FROM users_auth "
                                             "WHERE id_user = ?", (token_data,)).fetchone()
        if user_needs_confirmation is None:
            flash("User doesn't exist. Please contact PSCS admins.")
        if user_needs_confirmation['confirmed'] == 0:
            confirm_user(token_data)
            flash(f"User account for {user_needs_confirmation['name_user']} confirmed! Welcome to PSCS.")
        elif user_needs_confirmation['confirmed'] == 1:
            flash("User account has already been confirmed.")
    return redirect(url_for('index'))

def confirm_user(id_user):
    db = get_db()
    db.execute("UPDATE users_auth "
               "SET confirmed = 1, confirmed_datetime = CURRENT_TIMESTAMP "
               "WHERE id_user = ?", (id_user,))
    db.commit()
    return

@bp.route('login', methods=('GET', 'POST'))
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        user_login(username=username, password=password)
        return redirect(url_for('pscs.index'))
    return render_template('auth/login.html')


def user_login(username: str,
               password: str):
    """
    Checks whether the submitted password is valid for the username, and logs the user into the website.
    Parameters
    ----------
    username : str
        Username of the account
    password : str
        Password connected to the account

    Returns
    -------
    None
    """
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
        # check if user is admin
        is_admin = db.execute("SELECT id_user FROM admin_user WHERE id_user = ?", (user["id_user"],)).fetchone()
        if is_admin is not None:
            session["is_admin"] = True
        else:
            session["is_admin"] = False
    else:
        flash(error)
    return


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
