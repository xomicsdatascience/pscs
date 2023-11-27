from flask import session, g
from pscs.db import get_db
import uuid


def register_result(filename: str,
                    result_type: str,
                    title: str,
                    description: str):
    """
    SQL wrapper for registering a result into the results table.
    Parameters
    ----------
    filename : str
        Name of the file
    result_type : str
        Type of result (e.g., 'graph', 'table')
    title : str
        Title of the result, appropriate as a header.
    description : str
        Full description of the result.

    Returns
    -------
    None
    """
    db = get_db()
    result_row = 'start'
    while len(result_row) > 0:
        id_result = str(uuid.uuid4())
        result_row = db.execute('SELECT id_result FROM results WHERE id_result = ?', (id_result,)).fetchall()
    id_user = g.user['id_user']
    id_project = session['CURRENT_PROJECT']
    id_analysis = 0  # TODO
    # f = open(filename, 'rb')
    # file_hash = hashlib.sha3_256(f.read()).hexdigest()
    # f.close()

    db.execute('INSERT INTO results (id_result, id_project, id_analysis, file_path, result_type, title, description) VALUES (?,?,?,?,?,?,?)', (id_result, id_project, id_analysis, filename, result_type, title, description))
    db.commit()
    return