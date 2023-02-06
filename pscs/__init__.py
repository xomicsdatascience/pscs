__version__ = '0.0.1'

import os
from flask import Flask

def create_app(test_config=None) -> Flask:
    """
    Creates an app, configures it, and returns it.
    Parameters
    ----------
    test_config : mapping
        Flask mapping for testing a configuration.

    Returns
    -------
    Flask
        An instance of the Flask app.
    """
    app = Flask(__name__, instance_relative_config=True)
    # Recaptcha settings
    recaptcha_client_file = open(os.path.join(app.instance_path, 'recaptcha_client'), 'r')
    recaptcha_client_token = recaptcha_client_file.read().strip()
    recaptcha_client_file.close()
    recaptcha_server_file = open(os.path.join(app.instance_path, 'recaptcha_server'), 'r')
    recaptcha_server_token = recaptcha_server_file.read().strip()
    recaptcha_server_file.close()
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'pscs.sqlite'),
        RECAPTCHA_CLIENT=recaptcha_client_token,
        RECAPTCHA_SERVER=recaptcha_server_token)

    if test_config is None:
        # Load instance
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    # Make sure path exists
    os.makedirs(app.instance_path, exist_ok=True)

    from . import db
    db.init_app(app)

    from . import auth
    app.register_blueprint(auth.bp)

    from . import pscs
    app.register_blueprint(pscs.bp)
    app.add_url_rule('/', endpoint='index')

    return app
