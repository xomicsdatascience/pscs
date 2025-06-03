import sqlite3
import click
from flask import current_app, g
import json
import uuid
import os

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
        g.db.execute('PRAGMA foreign_keys = ON;')  # enforce foreign key relationships
    return g.db


def close_db(e=None):
    """Closes the currently-open db if one was opened."""
    db = g.pop('db', None)
    if db is not None:
        db.close()
    return


def init_db():
    # Raise error if already exists; this is to avoid overwriting it
    if os.path.exists(current_app.config["DATABASE"]):
        raise ValueError(f"Database already exists at {current_app.config['DATABASE']}")
    db = get_db()
    with current_app.open_resource('schema.sql') as f:
        db.executescript(f.read().decode('utf8'))
    # Populate table for university domains
    f = current_app.open_resource('universities.json')
    universities = json.load(f)
    for univ in universities:
        id_univ = get_unique_value_for_field(db, 'id_university', 'universities')
        db.execute('INSERT INTO universities (id_university, name, alpha_two_code, state_province, country, web_page)'
                   ' VALUES (?,?,?,?,?,?)',
                   (id_univ, univ['name'], univ['alpha_two_code'], univ['state-province'], univ['country'], univ['web_pages'][0]))
        for dom in univ['domains']:
            db.execute('INSERT INTO university_domains VALUES (?,?)', (id_univ, dom))
    for tag in current_app.config['INTERACTIVE_TAGS']:
        db.execute('INSERT INTO analysis_interactive_tags (interactive_tag) VALUES (?)', (tag,))
    db.commit()
    f.close()
    return


@click.command('init-db')
def init_db_command():
    """Clear existing data; create new tables."""
    init_db()
    click.echo('Database initialized.')
    click.echo(f'Location: {current_app.config["DATABASE"]}')
    return


def init_app(app):
    """Registers init db command"""
    app.teardown_appcontext(close_db)  # Call this function when shutting down
    app.cli.add_command(init_db_command)  # Add command to flask
    return


def get_unique_value_for_field(db, field: str, table: str) -> str:
    """
    Gets a unique uuid for the given field in the given table. This is intended for use with DBMs that don't have
    this as an built-in feature.
    Parameters
    ----------
    db
        Database object from which to pull, obtained via get_db
    field : str
        Field for which the unique value should be obtained
    table : str
        Table from which the field should be pulled

    Returns
    -------
    str
        Unique value for the field
    """
    row = 'temp'
    while len(row) != 0:
        tentative_id = str(uuid.uuid4())
        row = db.execute(f'SELECT {field} FROM {table} WHERE {field} = ?', (tentative_id,)).fetchall()
    return tentative_id


def delete_project(db: sqlite3.Connection, id_project: str):
    """
    Stages the project specified by the user ID and related content for deletion.
    Parameters
    ----------
    db : sqlite3.Connection
        DB from which to remove project.
    id_project : str
        ID of the project to delete.

    Returns
    -------
    None
    """
    db.execute('DELETE FROM projects WHERE id_project = ?', (id_project,))
    db.commit()
    return


def check_user_permission(permission_name: str,
                          permission_value: int,
                          id_project: str) -> bool:
    """
    Checks that the user has the appropriate permission for the specified project.
    Parameters
    ----------
    permission_name : str
        Name of the permission to check.
    permission_value
        Value that the permission should have to be accepted.
    id_project : str
        ID of the project to check.

    Returns
    -------
    bool
        Whether the current user has permission.
    """
    db = get_db()
    if permission_name in current_app.config["PROJECT_PERMISSIONS"]:  # prevent SQL injection
        role_info = db.execute(f'SELECT {permission_name} FROM projects_roles WHERE id_user = ? and id_project = ?',
                               (g.user['id_user'], id_project)).fetchone()
    else:
        return False
    if role_info is None:
        return False
    return role_info[permission_name] == permission_value


def check_analysis_published(id_analysis: str) -> bool:
    """Checks whether the analysis has been made publicly available."""
    db = get_db()
    # Check if default
    if len(db.execute("SELECT id_analysis FROM default_analysis WHERE id_analysis = ?", (id_analysis,)).fetchall()) > 0:
        return True
    analysis_info = db.execute("SELECT status FROM publications AS P "
                               "INNER JOIN publications_analysis AS PA on P.id_publication = PA.id_publication "
                               "WHERE PA.id_analysis = ?", (id_analysis,)).fetchall()
    if analysis_info is None:
        return False
    return any([a["status"] == "public" for a in analysis_info])


def check_analysis_is_default(id_analysis: str) -> bool:
    """Checks whether the analysis is one of the default analyses."""
    db = get_db()
    return db.execute("SELECT COUNT(*) FROM default_analysis WHERE id_analysis = ?", (id_analysis,)).fetchone()[0] > 0


def check_data_published(id_data: str) -> bool:
    """Checks whether the data has been made publicly available."""
    db = get_db()
    data_info = db.execute("SELECT status FROM publications AS P "
                           "INNER JOIN publications_data AS PD on P.id_publication = PD.id_publication "
                           "WHERE PD.id_data = ?", (id_data,)).fetchall()
    if data_info is None:
        return False
    return any([d["status"] == "public" for d in data_info])

def get_default_analysis() -> list:
    db = get_db()
    return db.execute("SELECT DA.id_analysis, A.analysis_name " 
                      "FROM default_analysis as DA INNER JOIN analysis AS A ON DA.id_analysis = A.id_analysis").fetchall()
