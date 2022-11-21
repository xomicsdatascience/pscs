import sqlite3
import click
from flask import current_app, g


def get_db():
    """
    Returns the database from currently-running app
    Returns
    -------
    sqlite3.database
        Current app's database.
    """
    if 'db' not in g:
        g.db = sqlite3.connect(
            current_app.config['DATABASE'],
            detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row

    return g.db


def close_db(e=None):
    """Closes the currently-open db if one was opened."""
    db = g.pop('db', None)
    if db is not None:
        db.close()
    return

def init_db():
    db = get_db()
    with current_app.open_resource('schema.sql') as f:
        db.executescript(f.read().decode('utf8'))
    return


@click.command('init-db')
def init_db_command():
    """Clear existing data; create new tables."""
    init_db()
    click.echo('Database initialized.')
    return


def init_app(app):
    """Registers init db command"""
    app.teardown_appcontext(close_db)  # Call this function when shutting down
    app.cli.add_command(init_db_command)  # Add command to flask