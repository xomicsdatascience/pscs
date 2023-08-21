"""This file contains code for resetting user passwords."""

from flask import current_app, flash, redirect
from pscs.messaging.mail import send_email
from pscs.db import get_db
from itsdangerous.url_safe import URLSafeTimedSerializer
from itsdangerous.exc import SignatureExpired, BadSignature


def send_reset_email(id_user):
    """Sends a password reset email to the user."""
    # Ge user email
    db = get_db()
    user_info = db.execute("SELECT email, name_user "
                           "FROM users_auth "
                           "WHERE id_user = ?", (id_user,)).fetchone()

    url_signer = URLSafeTimedSerializer(secret_key=current_app.config["SECRET_KEY"], salt="reset")
    token = url_signer.dumps(id_user)
    reset_url = "https://" + current_app.config["CURRENT_URL"] + f"/auth/reset/{token}"
    from pscs.templates.misc import password_reset
    password_formatted = password_reset.format(name_user=user_info["name_user"],
                                               url=reset_url)
    send_email([user_info["email"]], subject="PSCS - Password reset", body=password_formatted)
    return
