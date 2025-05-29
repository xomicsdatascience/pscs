import functools
import shutil
import sqlite3
import uuid

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
from werkzeug.utils import escape
import os
from os.path import join
import random
import pathlib

bp = Blueprint('auth', __name__, url_prefix='/auth')


@bp.route('/register', methods=('GET', 'POST'))
def register():
    """Create a new user."""
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        password_confirm = request.form['passwordConfirm']
        email = request.form['email']
        user_ip = request.remote_addr

        if current_app.config["RECAPTCHA_ENABLED"]:  # check that recaptcha is enabled
            recaptcha = request.form['g-recaptcha-response']
            recaptcha_valid, recaptcha_msg = validate_recaptcha(recaptcha, current_app.config['RECAPTCHA_SERVER'],user_ip)
            # If recaptcha is not valid, stop bothering the server
            if not recaptcha_valid:
                _ = recaptcha_msg  # Error from recaptcha is relatively useless; just tell user that there was a problem
                flash("Error with reCAPTCHA token.")
                return render_template("auth/register.html", recaptcha_enabled=current_app.config["RECAPTCHA_ENABLED"])
        db = get_db()
        try:
            return create_user(db, username, email, password, password_confirm, user_ip, request.form)
        except db.IntegrityError:
            error = f"User {username} is already registered."
        flash(error)
    return render_template("auth/register.html", recaptcha_enabled=current_app.config["RECAPTCHA_ENABLED"])


def create_user(db: sqlite3.Connection,
                username: str,
                email: str,
                password: str,
                password_confirm: str,
                user_ip: str,
                request_form: dict = None,
                is_temp_user: bool = False):
    # Generate uuid for user
    id_user = get_unique_value_for_field(db, 'id_user', 'users_auth')
    if request_form is None:
        if is_temp_user:
            request_form = {"noPHI": "on", "dataUse": "on"}
        else:
            raise ValueError("request_form must be defined for non-temp users.")

    is_registration_valid, err_msg = validate_registration(db, email, password, password_confirm, username, request_form, is_temp_user)
    if not is_registration_valid:
        flash(err_msg)
        return render_template("auth/register.html", recaptcha_enabled=current_app.config["RECAPTCHA_ENABLED"])
    else:
        if "name" not in request_form.keys():
            name = None
        else:
            name = request_form["name"]
        db.execute(
            "INSERT INTO users_auth (id_user, name_user, password, email, ip, confirmed_phi, confirmed_datause, name, is_temp_user) VALUES (?, ?, ?, ?, ?, ?, ?, ?,?)",
            (id_user, username, generate_password_hash(password), email, user_ip, request_form["dataUse"] == "on", request_form["noPHI"] == "on", name, is_temp_user))
        db.commit()

        # Send mail to user to verify.
        if current_app.config["REGISTRATION_REQUIRES_CONFIRMATION"] and not is_temp_user:
            send_user_confirmation_email(id_user, email, username)
            flash(f"Email confirmation sent to {email}; please click the verification link.")
        else:  # email confirmation disabled
            confirm_user(id_user)
        user_login(username, password)
        return redirect(url_for("pscs.index"))


def create_temp_user(db: sqlite3.Connection,
                     user_ip: str):
    username = make_temp_username(db)
    password = str(uuid.uuid4())
    email = get_unique_value_for_field(db, "email", "users_auth")
    _ = create_user(db, username = username, email = email, password = password, password_confirm=password, user_ip=user_ip, is_temp_user=True)
    load_logged_in_user()
    return


def make_temp_username(db: sqlite3.Connection):
    """Finds a unique username."""
    f = current_app.open_resource("static/wordlists/adjectives.txt", "r")
    adjectives = f.read().splitlines()
    f.close()
    f = current_app.open_resource("static/wordlists/nouns.txt", "r")
    nouns = f.read().splitlines()
    f.close()
    name_is_valid = False
    count = 0
    username = ""
    while not name_is_valid:
        i = random.randint(0, len(adjectives) - 1)
        j = random.randint(0, len(adjectives) - 1)
        while j == i:
            j = random.randint(0, len(adjectives) - 1)
        k = random.randint(0, len(nouns) - 1)
        adj1 = adjectives[i]
        adj2 = adjectives[j]
        noun = nouns[k]
        username = "".join(w.title() for w in [adj1, adj2, noun])
        name_is_valid, _ = validate_username(username, db)
        count += 1
        if count > 50:
            raise ValueError("Unable to find a unique username.")
    return username


def validate_registration(db: sqlite3.Connection,
                          email: str,
                          password: str,
                          password_confirm: str,
                          username: str,
                          request_form: dict,
                          is_temp_user: bool = False):
    """Checks whether the inputs are valid for the db; return True if valid, False if not. If not, also returns a
    Flash message."""
    uname_valid, uname_msg = validate_username(username, db)
    password_valid, password_msg = validate_password(password, password_confirm)
    email_valid, email_msg = validate_email(email, db)
    if not is_temp_user:
        noPHI_valid, noPHI_msg = validate_PHI(request_form)
        datause_valid, datause_msg = validate_datause(request_form)
    else:
        noPHI_valid, datause_valid = True, True

    error_msgs = []
    if not uname_valid:
        error_msgs.append(uname_msg)
    if not password_valid and not is_temp_user:
        error_msgs.append(password_msg)
    if not email_valid and not is_temp_user:
        error_msgs.append(email_msg)
    if not noPHI_valid and not is_temp_user:
        error_msgs.append(noPHI_msg)
    if not datause_valid and not is_temp_user:
        error_msgs.append(datause_msg)

    if len(error_msgs) > 0:
        return False, "; ".join(error_msgs)
    return True, ""


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
    create_default_project(db, id_user)
    db.commit()
    return


def create_default_project(db: sqlite3.Connection,
                           id_user: str,
                           source_id_project: str = "bb35db7c-7f9b-40bd-a274-4f722a2b739d"):
    """
    Creates the default project for new users to see PSCS functionality.

    Parameters
    ----------
    db : sqlite3.Connection
        Connection to the DB.
    id_user : str
        ID of the user for which to create the project.
    source_id_project : str, optional
        ID of the project.

    Returns
    -------
    None
    """
    # Copy project with id_project to new
    source_project_info = dict(db.execute("SELECT * FROM projects WHERE id_project = ?", (source_id_project,)).fetchone())
    source_project_info["id_user"] = id_user
    id_project = get_unique_value_for_field(db, "id_project", "projects")
    source_project_info["id_project"] = id_project

    if not _validate_sql_for_table(db, "projects", source_project_info):
        raise ValueError(f"Incorrect column name: {k} not valid.")
    insert_cols = ", ".join(list(source_project_info.keys()))
    placeholders = ", ".join(["?"] * len(source_project_info))
    sql_command = f"INSERT INTO projects ({insert_cols}) VALUES({placeholders})"
    db.execute(sql_command, tuple(source_project_info.values()))

    source_project_roles = dict(db.execute("SELECT * FROM projects_roles WHERE id_project = ?", (source_id_project,)).fetchone())
    source_project_roles["id_project"] = id_project
    source_project_roles["id_user"] = id_user
    insert_cols = ", ".join(list(source_project_roles.keys()))
    placeholders = ", ".join(["?"] * len(source_project_roles))
    sql_command = f"INSERT INTO projects_roles ({insert_cols}) VALUES({placeholders})"
    db.execute(sql_command, tuple(source_project_roles.values()))

    proj_dir = pathlib.Path(current_app.config['PROJECTS_DIRECTORY'].format(id_project=id_project))
    proj_dir.mkdir(exist_ok=True)
    data_dir = pathlib.Path(current_app.config["DATA_DIRECTORY"].format(id_project=id_project))
    data_dir.mkdir(exist_ok=True)
    results_dir = pathlib.Path(
        os.path.join(current_app.config['PROJECTS_DIRECTORY'], 'results').format(id_project=id_project))
    results_dir.mkdir(exist_ok=True)
    # Entry for project created; now copy data
    data_info = db.execute("SELECT * FROM data WHERE id_project = ?", (source_id_project,)).fetchall()
    data_info = [dict(d) for d in data_info]
    for d in data_info:
        d["id_data"] = get_unique_value_for_field(db, "id_data", "data")
        d["id_project"] = id_project
        new_path = join(current_app.config["DATA_DIRECTORY"].format(id_project=id_project), d["id_data"] + ".h5ad")
        shutil.copy(d["file_path"], new_path)
        d["file_path"] = new_path

        insert_cols = ", ".join(list(d.keys()))
        placeholders = ", ".join(["?"] * len(d))
        if not _validate_sql_for_table(db, "data", d):
            raise ValueError("Invalid data entry while creating default project.")
        sql_command = f"INSERT INTO data ({insert_cols}) VALUES({placeholders})"
        db.execute(sql_command, tuple(d.values()))

    # Copy analysis; import here to avoid circular imports
    from pscs.blueprints.projects import copy_pipeline
    analysis_info = db.execute("SELECT id_analysis FROM analysis WHERE id_project = ?", (source_id_project,)).fetchall()
    analysis_map = {}  # for results later
    for a in analysis_info:
        new_id = copy_pipeline(a["id_analysis"], id_project)
        analysis_map[a["id_analysis"]] = new_id

    # Copy results
    results_info = db.execute("SELECT * FROM results WHERE id_project = ?", (source_id_project,)).fetchall()
    amap_keys = list(analysis_map.keys())
    for rr in results_info:
        r = dict(rr)
        if r["id_result"] not in amap_keys:
            continue
        r["id_project"] = id_project
        previous_id = r["id_result"]
        r["id_result"] = get_unique_value_for_field(db, "id_result", "results")
        r["id_analysis"] = analysis_map[r["id_analysis"]]
        id_job = "id_job_sample"
        new_path = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=r["id_analysis"], id_job=id_job)
        previous_path_exp = os.path.basename(r["file_path"]).split(previous_id)[-1]  # get basename, remove id, get ext
        new_path = new_path + previous_path_exp
        r["file_path"] = new_path
        insert_cols = ", ".join(list(r.keys()))
        placeholders = ", ".join(["?"] * len(r))
        if not _validate_sql_for_table(db, "results", r):
            raise ValueError("Invalid result entry while creating default project.")
        sql_command = f"INSERT INTO results ({insert_cols}) VALUES({placeholders})"
        db.execute(sql_command, tuple(r.values()))

    db.commit()
    return


def _validate_sql_for_table(db: sqlite3.Connection,
                            table_name: str,
                            dict_to_validate: dict):
    """Checks whether the columns in `dict_to_validate` are in fact in the relevant table."""
    table_col_info = db.execute(f"PRAGMA table_info('{table_name}')").fetchall()
    col_names = set([c["name"] for c in table_col_info])
    for k in dict_to_validate.keys():
        if k not in col_names:
            return False
    return True


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
        'SELECT id_user, name_user, password, confirmed, is_temp_user FROM users_auth WHERE name_user = ?', (username,)
    ).fetchone()

    if user is None:
        error = 'Incorrect username.'
    elif not check_password_hash(user['password'], password):
        error = "Incorrect password."

    if error is None:
        # Deal with temporary user access
        tmp_user_id = None
        if "tmp_user_id" in session.keys():
            tmp_user_id = session["tmp_user_id"]

        session.clear()

        # Track temp ID for later integration into a real account
        if tmp_user_id is not None:
            session["tmp_user_id"] = tmp_user_id
        if user["is_temp_user"]:
            session["tmp_user_id"] = user["id_user"]

        session['id_user'] = user['id_user']
        session["confirmed"] = user["confirmed"]
        # check if user is admin
        is_admin = db.execute("SELECT id_user FROM admin_user WHERE id_user = ?", (user["id_user"],)).fetchone()
        if is_admin is not None:
            session["is_admin"] = True
        else:
            session["is_admin"] = False
        session.modified = True
        session.permanent = True
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
            "SELECT id_user, name_user, is_temp_user FROM users_auth WHERE id_user = ?", (user_id,)
        ).fetchone()


@bp.route('/logout')
def logout():
    if g.user is None:
        return redirect(url_for('index'))
    else:
        id_user = session["id_user"]
        is_temp = g.user["is_temp_user"]
        session.clear()
        if is_temp:
            try:
                delete_temp_user(id_user)
            except Exception as e:
                print(e)
    return redirect(url_for('index'))


def delete_temp_user(id_user):
    """Deletes a user flagged as temporary."""
    db = get_db()
    # Go specifically delete data
    data = db.execute("SELECT file_path FROM data WHERE id_user = ?", (id_user,)).fetchall()
    for d in data:
        os.remove(d["file_path"])
    db.execute("DELETE FROM users_auth WHERE id_user = ? AND is_temp_user = 1", (id_user,))
    db.commit()
    return


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


@bp.route("/updateName", methods=["POST"])
def update_name():
    """Updates the user's name"""
    if request.method == "POST":
        new_name = request.json["newName"]
        new_name = escape(new_name)
        if len(new_name) > current_app.config["MAX_NAME_LENGTH"] or len(new_name) < 2:
            flash("Name is too long; please contact PSCS admin.")
            return {"url": url_for("pscs.profile")}, 422
        else:
            db = get_db()
            db.execute("UPDATE users_auth SET name = ? WHERE id_user = ?", (new_name, g.user["id_user"]))
            db.commit()
            return {"url": url_for("pscs.profile")}, 200
    return {"url": url_for("pscs.profile")}, 405


@bp.route("/saveAffiliations", methods=["POST"])
def save_affiliations():
    """Updates the user's affiliations."""
    if request.method == "POST":
        affiliations = request.json["affiliations"]
        if len(affiliations) > current_app.config["MAX_AFFILIATIONS"]:
            flash("Too many affiliations; please reduce the number of affiliations and submit again.")
            return {"url": url_for("pscs.profile")}, 422
        if any([len(aff) > current_app.config["MAX_AFFILIATION_LENGTH"] for aff in affiliations]):
            flash("Affiliation names are too long; shorten")
            return {"url": url_for("pscs.profile")}, 422
        affiliations = [escape(aff) for aff in affiliations]
        # Remove previous affiliations; add new ones
        db = get_db()
        db.execute("DELETE FROM users_affiliation WHERE id_user = ?", (g.user["id_user"],))
        for idx_order, aff in enumerate(affiliations):
            db.execute("INSERT INTO users_affiliation (id_user, affiliation, affiliation_order) VALUES (?, ?, ?)",
                       (g.user["id_user"], aff, idx_order))
        db.commit()
        return {"url": url_for("pscs.profile")}, 200