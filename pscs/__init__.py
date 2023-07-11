__version__ = "0.0.49"

import os
from flask import Flask
import os
from os.path import join, dirname
import json
import shutil
import click


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
    env_dict = parse_env()
    static_directory = join(env_dict["INSTANCE_PATH"], "static")
    template_directory = join(env_dict["INSTANCE_PATH"], "templates")

    # Copy pscs package resources to instance
    this_dir = dirname(__file__)
    template_source = join(this_dir, "templates")
    static_source = join(this_dir, "static")
    print(f"Copying from {template_source} into {template_directory}")
    shutil.copytree(template_source, template_directory, dirs_exist_ok=True)
    shutil.copytree(static_source, static_directory, dirs_exist_ok=True)
    app = Flask(__name__, instance_path=env_dict["INSTANCE_PATH"],
                instance_relative_config=True, template_folder=template_directory,
                static_folder=static_directory)
    app.config.from_mapping(
        DATABASE=join(app.instance_path, 'pscs.sqlite'),
        **env_dict)
    if test_config is None:
        # Load instance
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    # Make sure path exists
    _makedir_until_format(app.instance_path)
    print(f"Instance path: {app.instance_path}")
    print(f"Making instance directories...")

    app.config["STATIC_DIRECTORY"] = static_directory
    _makedir_until_format(app.config["STATIC_DIRECTORY"])

    app.config['UPLOAD_FOLDER'] = join(app.instance_path, "upload", "{id_project}")
    _makedir_until_format(app.config['UPLOAD_FOLDER'])

    app.config["PROJECTS_DIRECTORY"] = join(app.instance_path, "projects", "{id_project}")
    _makedir_until_format(app.config["PROJECTS_DIRECTORY"])

    app.config["RESULTS_DIRECTORY"] = join(app.config["PROJECTS_DIRECTORY"], "results", "{id_analysis}")
    _makedir_until_format(app.config["RESULTS_DIRECTORY"])

    app.config["DELETION_DIRECTORY"] = join(app.instance_path, "deletion", "{id_project}")
    _makedir_until_format(app.config["DELETION_DIRECTORY"])

    app.add_url_rule('/upload/<name>', endpoint='pscs.download_file', build_only=True)

    app.cli.add_command(update_flask_files)
    from . import db
    db.init_app(app)

    from . import auth
    app.register_blueprint(auth.bp)

    from . import pscs
    app.register_blueprint(pscs.bp)

    from pscs.blueprints import homepage
    app.register_blueprint(homepage.bp)

    from pscs.blueprints import posts
    app.register_blueprint(posts.bp)

    from pscs.blueprints import projects
    app.register_blueprint(projects.bp)

    app.jinja_env.filters['basename'] = os.path.basename
    app.jinja_env.filters['dirname'] = os.path.dirname

    app.add_url_rule('/', endpoint='index')
    return app


@click.command('update-flask')
def update_flask_files():
    # Placeholder. The static files are updated when the app is called from the command line. If you want to set up
    # those files separately, move the relevant code here.
    return


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


def _makedir_until_format(path: str):
    """
    Makes the directories up until the directory contains a {format} tag.
    Parameters
    ----------
    path : str
        Path to make.

    Returns
    -------
    None
    """
    # Find {
    brace_idx = path.find("{")
    if brace_idx == -1:
        os.makedirs(path, exist_ok=True)
    else:
        os.makedirs(path[:brace_idx], exist_ok=True)
    return
