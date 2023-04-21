"""
This code is for handling email to/from the server. It is not intended to be used by users.
The code relies on having a Google Workspace account, setting up Google Cloud Services, enabling the Gmail API,
and getting the appropriate credentials + tokens.
Credentials/token paths can be changed.
"""
import os
import base64
import pickle
import json

# MIME messages
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

# Google OAuth2 + API
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError

token = "gmail_token.tok"
credentials_json = "gmail_json.json"

GOOGLE_API_ADDR = "https://www.googleapis.com/auth/gmail.send"


def _get_credentials_from_token(token_path: str):
    """Gets Gmail credentials from a stored token"""
    tok = open(token_path, "rb")
    tok_json = json.loads(pickle.load(tok))
    tok.close()
    credentials = Credentials.from_authorized_user_info(info=tok_json)
    return credentials


def _write_token(credentials,
                 token_path: str):
    """Overwrites the current token with the new one."""
    tok = open(token_path, "wb")
    pickle.dump(credentials.to_json(), tok)
    tok.close()
    return


def _get_credentials_from_secret(credentials_path: str):
    """
    USes the Google API with the app secret to get authentication token.
    Parameters
    ----------
    credentials_path : str
        Path to the app secrets JSON.
    Returns
    -------
    None
    """
    flow = InstalledAppFlow.from_client_secrets_file(credentials_path, [GOOGLE_API_ADDR])
    creds = flow.run_local_server(port=0)  # URL is printed to stdout; needs to be authorized
    return creds


def get_credentials():
    """Obtain credentials for app"""
    if os.path.exists(token):
        credentials =  _get_credentials_from_token(token)
        # Check validity; if expired and it can be refreshed,
        if not credentials.valid and credentials.expired and credentials.refresh_token:
            # Refresh credentials
            credentials.refresh(Request())
            _write_token(credentials, token)
    # No token; need to get from secret
    else:
        credentials = _get_credentials_from_secret(credentials_json)
        _write_token(credentials, token)
    return credentials


def send_email(recipients_to: list,
               subject: str,
               body: str,
               recipients_cc: list = None,
               recipients_bcc: list = None,
               ) -> bool:
    """
    Sends an email to the recipient.
    Parameters
    ----------
    recipients_to : list
        Recipients to list in the "To" field
    recipients_cc : list
        Recipients to list in the "cc" field
    recipients_bcc : list
        Recipients to list in the "bcc" field
    subject : str
        Subject line.
    body : str
        Content of the email.

    Returns
    -------
    bool
        Whether sending was successful.
    """
    # Create message
    message = MIMEMultipart()
    message["From"] = "PSCS <pscs@xods.org>"
    rec_to = ", ".join(recipients_to)
    message["To"] = rec_to
    if recipients_cc is not None:
        rec_cc = ", ".join(recipients_cc)
        message["Cc"] = rec_cc
    if recipients_bcc is not None:
        rec_bcc = ", ".join(recipients_bcc)
        message["Bcc"] = rec_bcc

    message["Subject"] = subject
    message.attach(MIMEText(body, "html"))

    # Get credentials
    credentials = get_credentials()
    service = build("gmail", "v1", credentials=credentials)
    message64 = {"raw": base64.urlsafe_b64encode(message.as_bytes()).decode()}
#               (service.users().messages().send(userId="me", body=create_message).execute())
    sent_info = (service.users().messages().send(userId="me", body=message64).execute())
    return "SENT" in sent_info["labelIds"]
