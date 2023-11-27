import functools
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app, jsonify
)

from werkzeug.security import check_password_hash, generate_password_hash
from itsdangerous.url_safe import URLSafeTimedSerializer
from pscs.db import get_db, get_unique_value_for_field
from pscs.authtools.validation.registration import validate_username, validate_password, validate_email, \
                                                   validate_recaptcha, send_user_confirmation_email, decode_token, \
                                                   validate_PHI, validate_datause
from pscs.authtools.validation.password_reset import send_reset_email
import math
import datetime
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
            noPHI_valid, noPHI_msg = validate_PHI(request.form)
            datause_valid, datause_msg = validate_datause(request.form)
            error = ''
            if not uname_valid:
                error = uname_msg + "; "
            if not password_valid:
                error += password_msg + "; "
            if not email_valid:
                error += email_msg + "; "
            if not noPHI_valid:
                error += noPHI_msg + "; "
            if not datause_valid:
                error += datause_msg + "; "
            if len(error) > 0:
                error = error[:-2]
                flash(error)
                return render_template("auth/register.html")
            else:
                db.execute(
                "INSERT INTO users_auth (id_user, name_user, password, email, ip, confirmed_phi, confirmed_datause) VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (id_user, username, generate_password_hash(password), email, user_ip, noPHI_valid, datause_valid))
                db.commit()

                # Send mail to user to verify.
                if current_app.config["REGISTRATION_REQUIRES_CONFIRMATION"]:
                    send_user_confirmation_email(id_user, email, username)
                    flash(f"Email confirmation sent to {email}; please click the verification link.")
                else:  # email confirmation disabled
                    confirm_user(id_user)
                user_login(username, password)
                return redirect(url_for("pscs.index"))
        except db.IntegrityError:
            error = f"User {username} is already registered."
        flash(error)
    return render_template("auth/register.html", recaptcha_enabled=current_app.config["RECAPTCHA_ENABLED"])


@bp.route('/confirmation/<token>', methods=["GET"])
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


@bp.route("/confirmation/resend/<token>", methods=["GET"])
def resend_confirmation(token):
    """Send confirmation email"""
    valid_token, token_data, invalid_reason = decode_token(token, "resend")
    if not valid_token:
        return redirect(url_for("pscs.indx"))
    db = get_db()
    user_info = db.execute("SELECT name_user, email "
                           "FROM users_auth "
                           "WHERE id_user = ?", (token_data,)).fetchone()
    if user_info is not None:
        send_user_confirmation_email(id_user=token_data, user_email=user_info["email"], name_user=user_info["name_user"])
    return redirect(url_for("pscs.index"))



@bp.route('/reset/<token>', methods=["GET", "POST"])
def password_reset(token):
    if request.method == "GET":
        valid_token, token_data, invalid_reason = decode_token(token, "reset", current_app.config["RESET_TIMEOUT_SECONDS"])
        if not valid_token and token != "debug":
            if len(invalid_reason) > 0:
                flash(invalid_reason)
                return redirect(url_for("pscs.index"))
        elif valid_token or token == "debug":
            # reset password
            return render_template("auth/password_reset.html", id_user=token_data)
    elif request.method == "POST":
        data = request.form
        password_valid, password_msg = validate_password(data["password"], data["passwordConfirm"])
        if not password_valid:
            flash(password_msg)
            return redirect(request.url)
        else:
            valid_token, token_data, invalid_reason = decode_token(token, "reset",
                                                                   current_app.config["RESET_TIMEOUT_SECONDS"])
            if not valid_token:
                flash(invalid_reason)
                return redirect(url_for("pscs.index"))
            elif valid_token:
                logout()
                db = get_db()
                db.execute("UPDATE users_auth "
                           "SET password = ? "
                           "WHERE id_user = ?", (generate_password_hash(data["password"]), token_data))
                db.commit()
                flash("Password updated.")
                return redirect(url_for("pscs.index"))


def confirm_user(id_user):
    db = get_db()
    db.execute("UPDATE users_auth "
               "SET confirmed = 1, confirmed_datetime = CURRENT_TIMESTAMP "
               "WHERE id_user = ?", (id_user,))
    session["confirmed"] = 1
    db.commit()
    return


@bp.route('login', methods=('GET', 'POST'))
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        err = user_login(username=username, password=password)
        if err is None:
            return redirect(url_for('pscs.index'))
        else:
            flash(err)
            redirect(url_for('auth.login'))
    return render_template('auth/login.html')


def user_login(username: str,
               password: str) -> str:
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
    str
        Error message. None if no error.
    """
    db = get_db()
    error = None
    user = db.execute(
        'SELECT id_user, name_user, password, confirmed FROM users_auth WHERE name_user = ?', (username,)
    ).fetchone()

    if user is None:
        error = 'Incorrect username.'
    elif not check_password_hash(user['password'], password):
        error = "Incorrect password."

    if error is None:
        session.clear()
        session['id_user'] = user['id_user']
        session["confirmed"] = user["confirmed"]
        # check if user is admin
        is_admin = db.execute("SELECT id_user FROM admin_user WHERE id_user = ?", (user["id_user"],)).fetchone()
        if is_admin is not None:
            session["is_admin"] = True
        else:
            session["is_admin"] = False
    return error


def check_password(id_user: str,
                   password: str) -> bool:
    """Checks the password against the DB"""
    db = get_db()
    user_pass = db.execute("SELECT password "
                           "FROM users_auth "
                           "WHERE id_user = ?", (id_user,))
    return check_password_hash(user_pass["password"], password)


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
        elif current_app.config["REGISTRATION_REQUIRES_CONFIRMATION"]:
            if "confirmed" not in session.keys() or session["confirmed"] != 1:
                # Check whether email has timed out
                db = get_db()
                id_user = g.user["id_user"]
                user_conf_info = db.execute("SELECT num_sent, last_sent "
                                            "FROM users_confirmation "
                                            "WHERE id_user = ?", (id_user,)).fetchone()
                is_okay_to_send = check_confirmation(user_conf_info)

                url_signer = URLSafeTimedSerializer(secret_key=current_app.config["SECRET_KEY"], salt="resend")
                token = url_signer.dumps(id_user)
                if is_okay_to_send:
                    resend_url = url_for("auth.resend_confirmation", token=token)
                    msg = f"Email verification is required. To send a new verification link, <a href='{current_app.config['CURRENT_URL']}{resend_url}'>click here</a>."
                else:
                    msg = "Email verification is required. Check your inbox for the confirmation email."
                flash(msg, category="link")
                return redirect(url_for('pscs.index'))
        return view(**kwargs)
    return wrapped_view


def is_logged_in():
    """Checks whether the user is logged in and confirmed (if confirmation is required)"""
    if g.user is None:
        return False
    elif current_app.config["REGISTRATION_REQUIRES_CONFIRMATION"]:
        if "confirmed" not in session.keys() or session["confirmed"] != 1:
            # Check whether email has timed out
            return False
    return True


def check_confirmation(user_conf_info: dict):
    """Checks whether it's okay to send another confirmation email."""
    if user_conf_info is None:
        return True
    cooldown = _get_confirmation_cooldown(user_conf_info["num_sent"])
    last_sent = user_conf_info["last_sent"]
    # Convert to epoch
    last_sent_epoch = last_sent.strftime("%s")
    last_sent_epoch = int(last_sent_epoch)
    # Convert current time to utc
    now_utc = int(datetime.datetime.now().astimezone(datetime.timezone.utc).strftime("%s"))
    time_since_sent = now_utc - last_sent_epoch
    return time_since_sent >= cooldown


def _get_confirmation_cooldown(num_sent: int) -> int:
    max_confirms = current_app.config["REGISTRATION_MAX_CONFIRMATION"]
    cooldown_time_minutes = current_app.config["REGISTRATION_CONFIRMATION_COOLDOWN"]
    if max_confirms > max_confirms:
        return math.inf  # too many attempts; label as spam
    else:
        return cooldown_time_minutes[num_sent]


def admin_required(view):
    @functools.wraps(view)
    def wrapped_view(**kwargs):
        if "is_admin" in session.keys() and not session["is_admin"]:
            return redirect(url_for("pscs.index"))
        return view(**kwargs)
    return wrapped_view


@bp.route('/reset/', methods=["GET", "POST"])
def reset_password_email():
    """Sends a password reset email to the user."""
    if request.method == "GET":
        return render_template("auth/password_forgot.html")
    if request.method == "POST":
        data = request.json
        id_user = None
        if "id_user" in data.keys():
            id_user = data["id_user"]
        elif "name_user" in data.keys():
            db = get_db()
            id_user = db.execute("SELECT id_user "
                                 "FROM users_auth "
                                 "WHERE name_user = ?", (data["name_user"],)).fetchone()["id_user"]
        if id_user is not None:
            send_reset_email(id_user)
            flash("The password reset information has been sent to your email.")
            return jsonify({"url": url_for("pscs.index")})
        else:
            flash("There was a problem retrieving the user email. Please contact PSCS admin.")
            return redirect(url_for("pscs.index"))
