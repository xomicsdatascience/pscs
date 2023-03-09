'''
This file contains the server-side validation for user registration.
'''
from pscs.db import get_db
from flask import current_app
import requests
import json
from itsdangerous.url_safe import URLSafeTimedSerializer
from itsdangerous.exc import SignatureExpired, BadSignature
from typing import Optional
import subprocess
min_password_length = 12


def validate_username(username: str, db=None) -> (bool, str):
    """
    Checks that the username exists, meets our requirements, and is not already in the database.
    Parameters
    ----------
    username : str
        Requested username
    db : sqlite3.Connection
        Optional. Connection to DB to check against. If not defined, doesn't check that username is unique.
    Returns
    -------
    bool
        Whether the username is valid
    str
        Message associated with the error, if any.
    """
    # NOTE: Other than the check for uniqueness, validation here only fails if user circumvented client-side
    # validation.
    msg = ''
    if len(username) == 0:
        msg = "Username is empty."
    elif not username[0].isalpha():
        msg = "First character must be alphabetical."
    elif db is not None:  # If no other errors, check against DB
        other_usernames_in_db = db.execute('SELECT name_user FROM users_auth WHERE name_user = ?', (username,)).fetchall()
        if len(other_usernames_in_db) > 0:
            msg = "Username already taken."
    return (len(msg) == 0), msg


def validate_email(email: str, db = None) -> (bool, str):
    """
    Checks that the email has the right form
    Parameters
    ----------
    email : str
        String containing the user's email address.
    db : sqlite3.Connection
        Optional. Connection to DB to check against. If not defined, email is not checked for uniqueness or for a valid
         domain.
    Returns
    -------
    bool
        Whether the email address is valid.
    str
        Message associated with email being invalid, if any.
    """
    msg = ''
    if len(email) == 0:
        msg = "Email address is empty."
    elif email.count('@') != 1:
        msg = "Email address format incorrect."
    elif db is not None:
        # Check for uniqueness
        db_emails = db.execute('SELECT email FROM users_auth WHERE email = ?', (email,)).fetchall()
        if len(db_emails) > 0:
            msg = "Email address is already used."
        # Validate domain
        domain = email.split('@')[-1]
        domain_split = domain.split('.')
        # Universities use subdomains for undergrad vs grad, departments, etc.
        # Those aren't part of our DB, so we're going to be less-specific until we get a match (if any)
        for idx in range(len(domain_split)):  # TODO: start from end of domain instead of the start; fewer sql pings
            test_domain = '.'.join(domain_split[idx:])
            db_domain = db.execute("SELECT university_domain FROM university_domains WHERE university_domain = ?", (test_domain,)).fetchall()
            if len(db_domain) == 0:
                continue
            break
        else:
            msg = "Email domain not associated with known university"
    return len(msg) == 0, msg


def validate_password(password: str, confirm_password: str) -> (bool, str):
    """
    Checks that the password and its confirmation are the same and that the password is adequate.
    Parameters
    ----------
    password : str
        Password sent by the user.
    confirm_password
        Password confirmation sent by the user.
    Returns
    -------
    bool
        Whether the password is valid.
    str
        Message associated with why the password is invalid, if any.
    """
    msg = ''
    if password != confirm_password:
        msg = "Password does not match confirmation."
    elif len(password) < min_password_length:
        msg = f"Password is too short ({min_password_length} min)"
    return len(msg) == 0, msg


def validate_recaptcha(recaptcha: str,
                       server_token: str,
                       user_ip: str = "0.0.0.0",
                       recaptcha_url: str = "https://www.google.com/recaptcha/api/siteverify") -> (bool, str):
    """
    Checks that the recaptcha token is valid.
    Parameters
    ----------
    recaptcha : str
        Client token obtained from the POST request sent by user.
    server_token : str
        Server secret token to use to confirm with recaptcha service.
    user_ip : str
        Requesting user's IP.
    recaptcha_url : str
        URL at which to validate the recaptcha

    Returns
    -------
    bool
        Whether the recaptcha is valid.
    str
        Message associated with why the recaptcha token is invalid, if any.
    """
    # Form POST content
    # secret, response, remoteip
    post_body = {"secret": server_token, "response": recaptcha, "remoteip": user_ip}
    response_txt = requests.post(recaptcha_url, post_body)
    response = json.loads(response_txt.text)
    if "error-codes" in response.keys():
        return response["success"], response["error-codes"]
    else:
        return response["success"], ""


def send_user_confirmation_email(id_user: str,
                                 user_email: str,
                                 name_user: str) -> bool:
    """
    Sends new user an email to confirm their address.
    Parameters
    ----------
    id_user : str
        User's ID.
    user_email : str
        User's email address.
    name_user : tuple
        Optional. Tuple of (First Name, Last Name) to use in email. If not used, those parts of the email are dropped.

    Returns
    -------
    bool
        Whether mail was sent successfully.
    """
    url_signer = URLSafeTimedSerializer(secret_key=current_app.config["SECRET_KEY"], salt="confirmation")
    token = url_signer.dumps(id_user)
    url = "https://" + current_app.config["CURRENT_URL"] + f"/auth/confirmation/{token}"
    message_contents = f"This email is to confirm the registration of user {name_user} with PSCS. Clicking the " \
                       f"<a href=\"{url}\">link</a> will confirm that you intended to register with PSCS. The link" \
                       f" will expire 12h from when it was issued.<br><br>"
    message_contents += "If you have received this message in error, no further action is needed on your part.<br><br>"
    message_contents += "- PSCS"
    return send_email(user_email, "PSCS confirmation", message_contents)


def send_email(email_address: str,
               subject: str,
               body: str) -> bool:
    """
    Send email body to address.
    Parameters
    ----------
    email_address : str
        Destination.
    subject : str
        Text to include in Subject line
    body : str
        Content of the message.

    Returns
    -------
    bool
        Whether mail was sent successfully.
    """
    msg = f"To: {email_address}\n"
    msg += f"Subject: {subject}\n"
    msg += "From: pscs@pscs.xods.org\n"
    msg += "MIME-Version: 1.0\n"
    msg += "Content-Type: text/html\n\n"  # multipart/alternative;
    msg += "<html>\n<body>\n"
    msg += body
    msg += "</body>\n</html>\n"
    sendmail_cmd = f"echo \"{msg}\" | sendmail {email_address}"
    return subprocess.run(sendmail_cmd, shell=True, capture_output=True)


def decode_token(token: str,
                 context: str,
                 max_age: Optional[int] = None) -> (bool, str):
    """
    Takes a token and verifies that it is valid with the specified context.
    Parameters
    ----------
    token : str
        Token to verify.
    context : str
        Context of the token. Used as salt with the most recent key.
    max_age : int, optional
        Maximum age in seconds that the token can be.
    Returns
    -------
    bool
        Whether the token is valid.
    object
        Data stored in the token.
    str
        Reason for invalid token, if any.
    """
    url_signer = URLSafeTimedSerializer(secret_key=current_app.config["SECRET_KEY"], salt=context)
    valid_token = False
    invalid_reason = None
    token_data = None
    try:
        token_data = url_signer.loads(token, max_age=max_age)
        valid_token = True
    except SignatureExpired:
        invalid_reason = "The token has expired."
    except BadSignature:
        invalid_reason = ""  # ignore the token
    return valid_token, token_data, invalid_reason
