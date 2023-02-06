'''
This file contains the server-side validation for user registration.
'''
from pscs.db import get_db
import requests
import json
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
    print(f"response: {response}")
    if "error-codes" in response.keys():
        return response["success"], response["error-codes"]
    else:
        return response["success"], ""


def send_user_confirmation_email(username: str) -> bool:
    """
    Sends new user an email to confirm their address.
    Parameters
    ----------
    username : str
        User's username.

    Returns
    -------
    bool
        Whether mail was sent successfully.
    """
    return True
