# This file contains the endpoints for the data uploader and associated code.
import os
import zipfile

from werkzeug.utils import secure_filename
from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app, send_from_directory,
    jsonify
)
import pandas as pd
import anndata as ad
from pscs.auth import login_required, is_logged_in
from pscs.pscs import delete_data
from pscs.db import get_unique_value_for_field, get_db, check_user_permission
from pscs.pscs import calc_hash_of_file
from os.path import join
from zipfile import ZipFile, ZIP_DEFLATED
import io

bp = Blueprint("uploader", __name__, url_prefix="/")


@bp.route("/<id_project>/upload/<uploader>", methods=["GET", "POST"])
@login_required
def uploader(id_project: str,
             uploader: str):
    """
    Returns the HTML page for the specified uploader
    Parameters
    ----------
    id_project : str
        ID of the project to which the data should be uploaded.
    uploader : str
        Name of the uploader to use.
    Returns
    -------
    str
        Rendered HTML
    """
    if request.method == "GET":
        if check_user_permission(permission_name="data_write",
                                 permission_value=1,
                                 id_project=id_project):
            # get the uploader
            url = url_for("projects.uploader.uploader", id_project=id_project, uploader=uploader)
            if uploader == "csv_h5ad":
                return render_template("pscs/uploader/csv_h5ad.html", form_url=url)
            elif uploader == "h5ad":
                return render_template("pscs/uploader/h5ad.html", form_url=url)
            else:
                flash("Uploader not found.")
                redirect(url_for("projects.upload", id_project=id_project))
        else:
            return redirect(url_for("pscs.index"))
    elif request.method == "POST":
        if check_user_permission(permission_name="data_write",
                                 permission_value=1,
                                 id_project=id_project):
            if uploader == "csv_h5ad":
                # try:
                process_csv_h5ad(request, id_project)
                flash("Data uploaded.")
                # except Exception as e:
                #     flash("There was a problem with the data upload.")
            elif uploader == "h5ad":
                process_h5ad(request, id_project)
                flash("Data uploaded.")
            return redirect(url_for("projects.project", id_project=id_project))
        else:
            return redirect(url_for("projects.project", id_project=id_project))


def extract_h5ad(req):
    """
    Given an HTML post request, returns
    Parameters
    ----------
    req : Flask.request
        HTML POST request sent by user.

    Returns
    -------
    flask.FileStorage
        Stream of the uploaded file; used for handling the incoming data.
    str
        Secure name of the uploaded file.
    """
    if req.form["file_name"] is None or len(req.form["file_name"]) == 0:
        fname = process_filename(req.files["h5ad"].filename)
    else:
        fname = process_filename(req.form["file_name"])
    return req.files["h5ad"].stream, fname


def process_h5ad(req, id_project: str):
    """
    Processes the request assuming that the uploaded file is an .h5ad file and registers it into the database.
    Parameters
    ----------
    req : Flask.request
        HTML POST request sent by user.
    id_project : str
        ID of the project to which the data should be uploaded.

    Returns
    -------
    None
    """
    h5ad_file_stream, h5ad_file_name = extract_h5ad(req)
    # No need to load data; just need to write it directly to disk.
    db = get_db()

    id_user = session["id_user"]
    id_data = get_unique_value_for_field(db, "id_data", "data")
    out_path = join(current_app.config["DATA_DIRECTORY"].format(id_project=id_project), id_data + ".h5ad")
    # Write out to datapath
    with open(out_path, "wb") as f:
        while True:
            chunk = h5ad_file_stream.read(16384)  # 16 kB read
            if not chunk:
                break
            f.write(chunk)
    file_hash = calc_hash_of_file(out_path)
    db.execute("INSERT INTO data (id_data, id_user, id_project, file_path, data_type, file_hash, file_name) "
               "VALUES (?,?,?,?,?,?,?)",
               (id_data, id_user, id_project, out_path, "h5ad", file_hash, h5ad_file_name))
    num_files = len(db.execute("SELECT id_data FROM data WHERE id_project = ?", (id_project,)).fetchall())
    db.execute(f'UPDATE projects SET num_files = {num_files} WHERE id_project = ?', (id_project,))

    # Don't need to make copy for original data; instead, place symlink in directory
    # (could do hardlink, but these might be on different filesystems in the future so meh)
    out_path_original_data = current_app.config["ORIGINAL_DATA_DIRECTORY"].format(id_project=id_project)
    os.makedirs(out_path_original_data, exist_ok=True)  # move to project creation
    id_data_original = get_unique_value_for_field(db, "id_data", "data_original")

    symlink_path = join(out_path_original_data, id_data_original)
    os.symlink(out_path, symlink_path)  # Create symlink pointing to out_path with name id_data_original.

    # Register - the fetcher is supposed to follow symlinks; check that it's doing that first
    db.execute(
        "INSERT INTO data_original (id_data, id_data_h5ad, id_user, id_project, file_path, data_type, file_name) "
        "VALUES (?,?,?,?,?,?,?)",
        (id_data_original, id_data, id_user, id_project, symlink_path, "h5ad", h5ad_file_name))
    db.commit()

    return


def process_filename(fname: str) -> str:
    """Processes uploaded filenames to be compatible with DB."""
    fname = secure_filename(fname)
    ext = os.path.splitext(fname)[1]
    if len(fname) > current_app.config["MAX_FILENAME_CHARS"]:
        fname = fname[:current_app.config["MAX_FILENAME_CHARS"] - len(ext)] + os.path.extsep + ext
    return fname[:current_app.config["MAX_FILENAME_CHARS"]]


def process_csv_h5ad(req, id_project: str):
    """
    Processes the uploaded data to be saved into PSCS.
    Parameters
    ----------
    req : Flask.request
        HTML POST request sent by user.
    id_project : str
    Returns
    -------

    """
    # First get data out of request
    data, uploaded_csv = extract_csv_h5ad(req)

    # Data ready to be written; get data ID
    db = get_db()
    # Get info: id_data, id_user, id_project, file_path, data_type, file_hash
    id_user = session["id_user"]
    if req.form["file_name"] is None or len(req.form["file_name"]) == 0:
        file_name = secure_filename(req.files["quantities"].filename)
    else:
        file_name = secure_filename(req.form["file_name"])
    # Remove file extension since it's been converted
    orig_file_name = file_name[:current_app.config["MAX_FILENAME_CHARS"]]
    last_ext = file_name.rfind(os.extsep)
    if last_ext != -1:
        file_name = file_name[:last_ext]
    file_name = file_name[:current_app.config["MAX_FILENAME_CHARS"]]
    data_type = "h5ad"
    id_data = get_unique_value_for_field(db, "id_data", "data")
    out_path = join(current_app.config["DATA_DIRECTORY"].format(id_project=id_project), id_data + ".h5ad")

    # AnnData (.h5ad) does not currently support __hash__
    # getting the hash otherwise requires us to convert the AnnData into bytes, which we don't want to support
    # Instead, save the file first, point hashlib at the file, get the hash
    data.write_h5ad(out_path, compression=current_app.config["DATA_COMPRESSION"])
    # Get hash
    file_hash = calc_hash_of_file(out_path)
    # Register data into db
    db.execute("INSERT INTO data (id_data, id_user, id_project, file_path, data_type, file_hash, file_name) "
               "VALUES (?,?,?,?,?,?,?)",
               (id_data, id_user, id_project, out_path, data_type, file_hash, file_name))
    num_files = len(db.execute("SELECT id_data FROM data WHERE id_project = ?", (id_project,)).fetchall())
    db.execute(f'UPDATE projects SET num_files = {num_files} WHERE id_project = ?', (id_project,))

    # Deal with uploaded data
    out_path_original_data = current_app.config["ORIGINAL_DATA_DIRECTORY"].format(id_project=id_project)
    os.makedirs(out_path_original_data, exist_ok=True)  # move to project creation
    id_data_original = get_unique_value_for_field(db, "id_data", "data_original")
    zip_filename = join(out_path_original_data, f"{id_data_original}.zip")
    with zipfile.ZipFile(zip_filename, 'w', ZIP_DEFLATED) as zipf:
        for k, v in uploaded_csv.items():
            csv_buffer = io.StringIO()
            v["data"].to_csv(csv_buffer, index=True)
            zipf.writestr(v["name"], csv_buffer.getvalue())
    # Register
    db.execute("INSERT INTO data_original (id_data, id_data_h5ad, id_user, id_project, file_path, data_type, file_name) "
               "VALUES (?,?,?,?,?,?,?)", (id_data_original, id_data, id_user, id_project, zip_filename, "csv", orig_file_name + ".zip"))
    db.commit()
    return


def extract_csv_h5ad(req) -> (ad.AnnData, dict):
    """
    Accepts three files: "quantities", "obs", and "var", and combines them into a single .h5ad file.
    Parameters
    ----------
    req
        Request to be processed.
    Returns
    -------
    ad.AnnData
        AnnData created from the uploaded data.
    dict
        pd.DataFrame
    """
    # Check what data is available
    quant = req.files["quantities"]
    do_quant_transpose = req.form.get("quant_transpose", None) is not None

    obs = req.files["obs_input"]
    do_obs_transpose = req.form.get("obs_transpose", None) is not None

    var = req.files["var_input"]
    do_var_transpose = req.form.get("var_transpose", None) is not None

    # Load data
    uploaded_csv = {}
    data_csv = pd.read_csv(quant.stream, index_col=0)   # add variable to first_column_names for data without idx
    uploaded_csv["quantities"] = {"name": secure_filename(quant.filename), "data": data_csv}
    data = ad.AnnData(X=data_csv)
    if do_quant_transpose:
        data = data.T

    if obs.filename != "":
        obs_data = pd.read_csv(obs.stream, index_col=0)  # index_col to false for data without idx
        if do_obs_transpose:
            obs_data = obs_data.T
        data.obs = obs_data
        uploaded_csv["obs"] = {"name": secure_filename(obs.filename), "data": obs_data}
    if var.filename != "":
        var_data = pd.read_csv(var.stream, index_col=0)  # index_col to false for data without idx
        if do_var_transpose:
            var_data = var_data.T
        data.var = var_data
        uploaded_csv["var"] = {"name": secure_filename(var.filename), "data": var_data}

    # ensure that filenames are unique
    fnames = set([v["name"] for v in uploaded_csv.values()])
    if len(fnames) != len(uploaded_csv):
        for k, v in uploaded_csv:
            fname = os.path.splitext(v["name"])
            v["name"] = fname[0] + f"_{k}" + fname[1]

    data.obs_names_make_unique()
    return data, uploaded_csv


@bp.route("/<id_project>/upload/<uploader>/script", methods=["GET"])
@login_required
def uploader_script(id_project: str,
                    uploader: str):
    """
    Returns the HTML page for the specified uploader
    Parameters
    ----------
    id_project : str
        ID of the project to which the data should be uploaded.
    uploader : str
        Name of the uploader to use.
    Returns
    -------
    str
        Rendered HTML
    """
    if check_user_permission(permission_name="data_write",
                             permission_value=1,
                             id_project=id_project):
        if uploader == "csv_h5ad":
            scripts = dict()
            scripts["js"] = ["/static/js/uploader/csv_h5ad.js", "/static/js/external/papaparse/papaparse.min.js"]
            return jsonify(scripts)
        elif uploader == "h5ad":
            scripts = dict()
            scripts["js"] = ["/static/js/uploader/h5ad.js"]
            return jsonify(scripts)
    return
