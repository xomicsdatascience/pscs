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
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'scpp.sqlite'))

    if test_config is None:
        # Load instance
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    # Make sure path exists
    os.makedirs(app.instance_path, exist_ok=True)

    # Say hello!
    @app.route('/hello')
    def hello():
        return "Hello!"

    from . import db
    db.init_app(app)

    from . import auth
    app.register_blueprint(auth.bp)

    from . import scpp
    app.register_blueprint(scpp.bp)
    app.add_url_rule('/', endpoint='index')

    return app