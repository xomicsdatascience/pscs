# This file contains the endpoints for the data uploader and associated code.
import os
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
            if uploader == "csv_h5ad":
                url = url_for("projects.uploader.uploader", id_project=id_project, uploader=uploader)
                return render_template("pscs/uploader/csv_h5ad.html", form_url=url)
        else:
            return redirect(url_for("pscs.index"))
    elif request.method == "POST":
        if uploader == "csv_h5ad":
            # try:
            process_csv_h5ad(request, id_project)
            flash("Data uploaded.")
            # except Exception as e:
            #     flash("There was a problem with the data upload.")
        return redirect(url_for("projects.project", id_project=id_project))


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
    data = extract_csv_h5ad(req)

    # Data ready to be written; get data ID
    db = get_db()
    # Get info: id_data, id_user, id_project, file_path, data_type, file_hash
    id_user = session["id_user"]
    if req.form["file_name"] is None or len(req.form["file_name"]) == 0:
        file_name = secure_filename(req.files["quantities"].filename)
    else:
        file_name = secure_filename(req.form["file_name"])
    # Remove file extension since it's been converted
    last_ext = file_name.rfind(os.extsep)
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
    db.commit()
    return


def extract_csv_h5ad(req):
    """
    Accepts three files: "quantities", "obs", and "var", and combines them into a single .h5ad file.
    Parameters
    ----------
    req
        Request to be processed.
    Returns
    -------

    """
    # Check what data is available
    quant = req.files["quantities"]
    do_quant_transpose = req.form.get("quant_transpose", None) is not None

    obs = req.files["obs_input"]
    do_obs_transpose = req.form.get("obs_transpose", None) is not None

    var = req.files["var_input"]
    do_var_transpose = req.form.get("var_transpose", None) is not None

    # Load data
    # data = ad.read_csv(quant.stream, first_column_names=True)
    data = pd.read_csv(quant.stream, index_col=0)   # add variable to first_column_names for data without idx
    data = ad.AnnData(X=data)
    if do_quant_transpose:
        data = data.T

    if obs.filename != "":
        obs_data = pd.read_csv(obs.stream, index_col=0)  # index_col to false for data without idx
        if do_obs_transpose:
            obs_data = obs_data.T
        data.obs = obs_data
        del obs_data  # clear memory asap
    if var.filename != "":
        var_data = pd.read_csv(var.stream, index_col=0)  # index_col to false for data without idx
        if do_var_transpose:
            var_data = var_data.T
        data.var = var_data
        del var_data  # clear memory asap
    data.obs_names_make_unique()
    return data


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
    return
