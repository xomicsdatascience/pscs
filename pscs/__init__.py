__version__ = "0.0.6"

import os
from flask import Flask
import os
import json
from waitress import serve


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
    app = Flask(__name__, instance_relative_config=True, template_folder="templates")
    env_dict = parse_env()
    app.config.from_mapping(
        DATABASE=os.path.join(app.instance_path, 'pscs.sqlite'),
        **env_dict)
    if test_config is None:
        # Load instance
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    # Make sure path exists
    os.makedirs(app.instance_path, exist_ok=True)
    print(f"Instance path: {app.instance_path}")

    from . import db
    db.init_app(app)

    from . import auth
    app.register_blueprint(auth.bp)

    from . import pscs
    app.register_blueprint(pscs.bp)

    from pscs.blueprints import homepage
    app.register_blueprint(homepage.bp)

    app.add_url_rule('/', endpoint='index')
    return app

def parse_env(env_file: str = '.env') -> dict:
    """
    Parses the environment file and returns the result as a dict.
    Parameters
    ----------
    env_file : str
        Path to the environment

    Returns
    -------
    dict
        Dictionary keyed with the environment variable name.
    """
    f = open(env_file)
    env_dict = dict()
    line = f.readline()
    line_count = 0  # on the off chance that the problematic line contains the secret key, don't print contents
    while len(line) > 0:
        line_count += 1
        idx_equal = line.index('=')
        if idx_equal == -1:
            raise ValueError(f"Environment file not correctly configured; see line {line_count} in .env file.")
        key = line[:idx_equal]
        env_dict[key] = json.loads(line[idx_equal+1:])
        line = f.readline()
    f.close()
    return env_dict
