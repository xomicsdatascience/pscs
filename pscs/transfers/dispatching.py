# This file contains code for dispatching processing pipelines to different computational resources.
import sqlite3
import subprocess

import anndata as ad
import pandas as pd
import numpy as np

from pscs.db import get_db, get_unique_value_for_field
from flask import current_app
from os.path import join, basename
import json
import tempfile
import os
import time
import threading


def dispatch(pipeline_json: str,
             id_user: str,
             file_info: dict,
             id_project: str,
             id_analysis: str,
             mem_requirement_mb: int = None,
             resource: str = 'osp') -> str:
    """
    Collects the relevant files and sends them to the relevant computational resource.
    Parameters
    ----------
    pipeline_json : str
        Path to the JSON file containing the node information.
    id_user : str
        User id of the submitter.
    file_info : dict
        Dict of {node_id: {"id" : id_data, "file_path": file_path}} containing the input files to use for the input nodes
         specified in pipeline_json.
    id_project : str
        ID of the project for which the analysis is being run.
    id_analysis : str
        ID of the analysis being run.
    mem_requirement_mb : int
        Number of MB of memory required.
    resource : str
        Name of the resource that should be used.
    Returns
    -------
    str
        PSCS job ID.
    """
    file_ids = dict()
    file_paths = dict()
    for k, v in file_info.items():
        file_ids[k] = file_info[k]["id"]
        file_paths[k] = file_info[k]["file_path"]

    db = get_db()
    pscs_job_id = get_unique_value_for_field(db, "id_job", "submitted_jobs")

    # Make destination
    remote_proj_dir = get_remote_project_dir(id_project=id_project, id_analysis=id_analysis, resource=resource)
    remote_mkdir(remote_dir=remote_proj_dir, resource=resource)
    remote_results = join(remote_proj_dir, 'results', pscs_job_id)
    remote_mkdir(remote_dir=remote_results, resource=resource)

    tags = db.execute("SELECT interactive_tag FROM analysis_interactive_tags").fetchall()
    tags = [tag["interactive_tag"] for tag in tags]
    for tag in tags:
        remote_tag_dir = join(remote_results, tag)
        remote_mkdir(remote_dir=remote_tag_dir, resource=resource)

    # Transfer files
    transfer_file(pipeline_json, remote_dir=remote_proj_dir, resource=resource)
    transfer_files(local_files=list(file_paths.values()), remote_dir=remote_proj_dir, resource=resource)

    # Create and transfer dictionary for new filepaths on remote machine
    remote_paths = remap_filepaths(file_paths, remote_dir=remote_proj_dir)
    remap_filename = save_filepath_remap(remote_paths)
    transfer_file(remap_filename, remote_dir=remote_proj_dir, resource=resource)

    # Submission script

    htcondor_script = generate_htcondor_submission(pipeline_json=pipeline_json,
                                                   remapped_input_json=remap_filename,
                                                   remote_project_dir=remote_proj_dir,
                                                   output_dir=remote_results,
                                                   id_user=id_user,
                                                   pscs_job_id=pscs_job_id,
                                                   mem_requirement_mb=mem_requirement_mb)
    transfer_file(htcondor_script, remote_proj_dir, resource=resource)
    cmd = f"condor_submit {os.path.join(remote_proj_dir, basename(htcondor_script))}"
    server_response = remote_cmd(cmd)

    # Log into db
    resource_job_id = get_job_id(server_response, resource)
    db.execute("INSERT INTO submitted_jobs "
               "(id_job, submitted_resource, resource_job_id, id_user, id_project, id_analysis, server_response, remote_results_directory) "
               "VALUES (?,?,?,?,?,?,?,?)", (pscs_job_id, resource, resource_job_id, id_user, id_project, id_analysis, server_response, remote_results))
    # Store submitted data in db
    for id_node, id_data in file_ids.items():
        db.execute("INSERT INTO submitted_data "
                   "(id_job, id_data, node_name) "
                   "VALUES (?, ?, ?)", (pscs_job_id, id_data, id_node))
    db.commit()
    return pscs_job_id


def local_dispatch(pipeline_json: str,
                   id_user: str,
                   file_info: dict,
                   id_project: str,
                   id_analysis: str,
                   resource: str = 'local'):
    if resource != "local":
        raise ValueError("local_dispatch can only be used with local resource")
    file_ids = dict()
    file_paths = dict()
    for k, v in file_info.items():
        file_ids[k] = file_info[k]["id"]
        file_paths[k] = file_info[k]["file_path"]
    db = get_db()
    pscs_job_id = get_unique_value_for_field(db, "id_job", "submitted_jobs")

    # Make local directories for files, results
    dirs_to_make = []
    proj_dir = current_app.config["PROJECTS_DIRECTORY"].format(id_project=id_project)
    dirs_to_make.append(proj_dir)
    results_dir = current_app.config["RESULTS_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis, id_job=pscs_job_id)
    dirs_to_make.append(results_dir)
    logs_dir = current_app.config["LOG_DIRECTORY"].format(id_project=id_project, id_analysis=id_analysis, id_job=pscs_job_id)
    print(f"Logs dir: {logs_dir}")
    dirs_to_make.append(logs_dir)
    tags = db.execute("SELECT interactive_tag FROM analysis_interactive_tags").fetchall()
    tags = [tag["interactive_tag"] for tag in tags]
    for tag in tags:
        tag_dir = join(results_dir, tag)
        dirs_to_make.append(tag_dir)

    make_local_directories(dirs_to_make)
    # file info
    f = open(f"{results_dir}/pscs-filepaths-{pscs_job_id}.json", 'w')
    json.dump(file_paths, f)
    f.close()
    print(f"\n\nFile info: {proj_dir}/pscs-filepaths-{pscs_job_id}.json\n\n")

    db.execute("INSERT INTO submitted_jobs "
               "(id_job, submitted_resource, resource_job_id, id_user, id_project, id_analysis, server_response, remote_results_directory) "
               "VALUES (?,?,?,?,?,?,?,?)", (pscs_job_id, resource, pscs_job_id, id_user, id_project, id_analysis, "", results_dir))
    db.execute("INSERT INTO local_jobs "
               "(id_job) VALUES (?)", (pscs_job_id,))
    db.commit()
    command = ["python", "pscs/run_local_pipeline.py", current_app.config['DATABASE'], current_app.config["REMOTE_COMPUTING"][resource]["URL"]]
    stdout_path = f"{logs_dir}/pscs-stdout-{pscs_job_id}.log"
    stderr_path = f"{logs_dir}/pscs-stderr-{pscs_job_id}.log"
    stdout_file = open(stdout_path, 'w')
    stderr_file = open(stderr_path, 'w')
    command = ["python", "pscs/run_local_pipeline.py", current_app.config['DATABASE'],
               current_app.config["REMOTE_COMPUTING"][resource]["URL"], stdout_path, stderr_path]
    process = subprocess.Popen(command, stdout=stdout_file, stderr=stderr_file)

    def cleanup():
        process.wait()
        stdout_file.close()
        stderr_file.close()

    # Start the cleanup in a background thread.
    cleanup_thread = threading.Thread(target=cleanup)
    cleanup_thread.start()

    return pscs_job_id

def make_local_directories(local_dirs: list):
    if isinstance(local_dirs, str):
        local_dirs = [local_dirs]
    for local_dir in local_dirs:
        os.makedirs(local_dir, exist_ok=True)
    return


def determine_resource(file_info: dict = None,
                       pipeline_json: str = None):
    remote_resources = current_app.config["REMOTE_COMPUTING"]
    # convert remote res dict to list; add name
    remote_resource_list = []
    for resource_name, resource_info in remote_resources.items():
        remote_resource_list.append(resource_info.copy())
        remote_resource_list[-1]["name"] = resource_name[0]
    sorted_resources = sorted(remote_resource_list, key=lambda x: x["preference"])

    # Estimate job memory requirements
    mem_requirement = estimate_job_memory_requirement(file_info=file_info,
                                                      pipeline_json=pipeline_json)
    mem_requirement = mem_requirement // (2**20)  # convert bytes to MB

    # Check which resources are eligible
    eligible_resources = []
    for res in sorted_resources:
        if res["max_memory_mb"] >= mem_requirement or res["max_memory_mb"] < 0:
            eligible_resources.append(res)

    # Check duration
    duration = np.inf
    selected_resource = None
    for res in eligible_resources:
        est_time = estimate_queue_length(res["name"])
        if est_time < duration:
            duration = est_time
            selected_resource = res["name"]
    if selected_resource is not None:
        return selected_resource
    else:
        # Something has gone wrong; submit to osp anyway
        return "osp", min(mem_requirement, 20*(2**1024))  # 20GB limit


def estimate_job_memory_requirement(file_info: dict, pipeline_json: str):
    """Estimates the memory requirement for a particular job given the input data and pipeline configuration."""
    max_mem = estimate_file_max_memory_requirements([f["file_path"] for f in file_info.values()])
    mult = estimate_pipeline_memory_multiplier(pipeline_json)
    return max_mem * mult


def estimate_file_max_memory_requirements(fpaths: list[str]):
    """Given a set of files, find the maximum estimated memory requirement across them."""
    # The value is multiplied by the number of paths; taking the max is conservative
    max_mem = None
    for fpath in fpaths:
        if fpath.endswith(".h5ad"):
            mem = estimate_h5ad_memory_reqs(fpath)
        else:
            mem = os.path.getsize(fpath)  # assumes no compression
        if max_mem is None:
            max_mem = mem
        else:
            max_mem = max(max_mem, mem)
    if max_mem is None:
        print("Unable to estimate memory requirements.")
        return 5000
    return max_mem


def estimate_pipeline_memory_multiplier(pipeline_json: str):
    """
    Estimates the memory multiplier for a pipeline. This estimate is aggressive and assumes that pipeline branches
    are never cleared. It also does not take into account that different paths may use different input files.
    Parameters
    ----------
    pipeline_json : str
        Path to the pipeline JSON file.
    Returns
    -------
    float
        Multiplier for memory requirements.
    """
    # Get number of splitting paths
    # NOTE: A better way would be to execute all nodes down to an output node, clear the object, then
    # repeat for all paths. The multiplier would be then be the maximum number of inputs of any node instead.
    mult = 0
    with open(pipeline_json, "r") as f:
        pjson = json.load(f)
        # Get number of inputs first
        for n in pjson["nodes"]:
            if n["num_inputs"] == 0:
                mult += 1
        for n in pjson["nodes"]:
            if n["num_inputs"] == 0:
                continue  # skip output nodes
            mult *= n["num_outputs"]
    return mult

def estimate_h5ad_memory_reqs(h5ad_filepath):
    """Estimates the memory requirements for loading an h5ad file into memory."""
    adata = ad.read_h5ad(h5ad_filepath, backed="r")
    adata_groups = ["X", "obs", "var", "obsm", "varm", "uns", "obsp", "varp"]
    adata_groups.extend(["_" + k for k in adata_groups])
    adata_dict = adata.__dict__
    mem_req = 0  # in bytes
    for group in adata_groups:
        if group in adata_dict:
            gdata = adata_dict[group]
            try:
                if group == "X":
                    mem_req += np.prod(gdata.shape) * gdata.dtype.itemsize
                # Check if df
                elif isinstance(gdata, pd.DataFrame):
                    mem_req += gdata.memory_usage(deep=False).sum()
                elif isinstance(gdata, ad._core.aligned_mapping.AxisArrays):
                    for k in gdata.keys():
                        mem_req += np.prod(gdata[k].shape) * gdata[k].dtype.itemsize
            except AttributeError:
                continue  # errors don't matter much here; it'll just result in a worse estimate
    return mem_req


def estimate_queue_length(resource):
    if resource == "local":
        return local_queue_length()
    elif resource == "osp":
        return 15*60  # roughly 15-minute jobs
    else:
        return np.inf


def local_queue_length():
    """Estimates the length of the local queue and returns it in seconds."""
    # We would need a lot of performance info for each node in a pipeline to get decent estimates
    # We'll instead get the number of jobs and assume ~60s per job.
    db = get_db()
    num_jobs = db.execute("SELECT COUNT(*) AS num_jobs FROM submitted_jobs "
                          "WHERE submitted_resource = 'local' AND is_complete = 0").fetchone()["num_jobs"]
    if num_jobs is None:
        return 0
    return 60*num_jobs


def is_local_queue_too_long():
    db = get_db()
    # get local jobs that are not completed
    jobs_longer_than_wait = db.execute("SELECT * "
                                 "FROM submitted_jobs "
                                 "WHERE (strftime('%s', 'now') - strftime('%s', date_submitted)) > ? "
                                 "AND submitted_resource = 'local' AND is_complete = 0",
                                 (current_app.config["LOCAL_MAX_WAIT_SECONDS"],)).fetchall()
    return len(jobs_longer_than_wait) > 0


def get_job_id(response: str,
               resource: str) -> str:
    """
    Extracts the job ID on the remote resource.
    Parameters
    ----------
    response : str
        Response from the server containing the job ID.
    resource : str
        Resource to which the job was submitted
    Returns
    -------
    str
        Job ID.
    """
    job_id = ""
    if resource == "osp":
        job_id = get_osp_job_id(response)
    return job_id


def get_osp_job_id(response: str) -> str:
    """Extracts the job ID from OSP's response"""
    # Format: "Submitting job(s).\nN jo(s) submitted to cluster JOB_ID
    split_response = response.split(' ')
    return split_response[-1].split(".")[0]


def save_filepath_remap(remote_paths: dict) -> str:
    """
    Saves the dictionary to a temporary file.
    Parameters
    ----------
    remote_paths : dict
        As returned from remap_filepaths

    Returns
    -------
    str
        Name of the temporary file on the local machine.
    """
    fd, filename = tempfile.mkstemp(prefix='pscs-filepaths-', suffix='.json')
    remap_file = os.fdopen(fd, 'w')
    json.dump(remote_paths, remap_file)
    remap_file.close()
    return filename


def remap_filepaths(file_paths: dict,
                    remote_dir: str) -> dict:
    """
    Converts the input file paths to paths on the remote resource.
    Parameters
    ----------
    file_paths : dict
        Dict of {node_id: file_path} containing the input files to use for the input nodes specified in pipeline_json.
    remote_dir : str
        Path on the remote resource.

    Returns
    -------
    dict
        Dict of {node_id: file_path} with the file_path updated to reflect the path on the remote resource.
    """
    remote_paths = {}
    for node_id, file_path in file_paths.items():
        # remote_paths[node_id] = join(remote_dir, basename(file_path))
        remote_paths[node_id] = basename(file_path)
    return remote_paths


def get_address(resource: str) -> str:
    """
    Returns the [user]@[url] for a given resource.
    Parameters
    ----------
    resource : str
        Resource for which to get the address.

    Returns
    -------
    str
        Address to use for resource
    """
    remote_user = current_app.config["REMOTE_COMPUTING"][resource]["USER"]
    remote_url = current_app.config["REMOTE_COMPUTING"][resource]["URL"]
    return f"{remote_user}@{remote_url}"


def remote_cmd(cmd: str,
               resource: str = 'osp'):
    """
    Executes the command on the remote resource.
    Parameters
    ----------
    cmd : str
        Command to run
    resource : str
        String specifying the resource on which to run the command. Default: 'osp'

    Returns
    -------
    None
    """
    if '$' in cmd:  # disallow variable expansion in the command
        return
    remote_address = get_address(resource)
    ssh_cmd = f"ssh {remote_address}"
    output = subprocess.check_output(ssh_cmd.split(), input=bytes(cmd, 'utf-8'))
    return output.decode('utf-8')


def get_remote_project_dir(id_project: str,
                           id_analysis: str,
                           resource: str = 'osp') -> str:
    """

    Parameters
    ----------
    resource : str
        Remote resource to use.
    id_project : str
        ID of project for which the analysis is being run.
    id_analysis : str
        ID of the analysis that is being run.

    Returns
    -------
    str
        Path on remote resource.
    """
    remote_outdir = join(current_app.config["REMOTE_COMPUTING"][resource]["OUTDIR"], id_project, id_analysis)
    return remote_outdir


def remote_mkdir(remote_dir: str,
                 resource: str = 'osp'):
    """
    Creates the specified directory in the remote resource.
    Parameters
    ----------
    remote_dir : str
        Directory to create on the remote resource.
    resource : str
        Remote resource to use.

    Returns
    -------
    None
    """
    remote_cmd(f'mkdir -p {remote_dir}', resource=resource)
    return


def transfer_file(local_file: str,
                  remote_dir: str,
                  resource: str = 'osp'):
    """
    Transfers file from local system to remote resource.
    Parameters
    ----------
    local_file : str
        Path to file to transfer.
    remote_dir : str
        Destination directory.
    resource : str
        Resource to which to transfer file.
    Returns
    -------
    None
    """
    remote_address = get_address(resource)
    cmd = f'rsync {local_file} {remote_address}:{remote_dir}'
    _ = subprocess.check_output(cmd.split())
    return


def transfer_files(local_files: list,
                   remote_dir: str,
                   resource: str = 'osp'):
    """
    Transfers multiple files from local system to remote resource.
    Parameters
    ----------
    local_files : list
        List of file paths to transfer
    remote_dir : str
        Directory at destination.
    resource : str
        Resource to which to transfer files.

    Returns
    -------
    None
    """
    for file in local_files:
        transfer_file(file, remote_dir, resource)
    return


def generate_htcondor_submission(pipeline_json: str,
                                 remapped_input_json: str,
                                 remote_project_dir: str,
                                 output_dir: str,
                                 id_user: str,
                                 pscs_job_id: str,
                                 mem_requirement_mb: int = 5000) -> str:
    """
    Creates a .submit script for resources using HTCondor.
    Parameters
    ----------
    pipeline_json : str
        Path to the JSON file containing the node information.
    remapped_input_json : str
        Path to JSON keyed with {node_id: file_path} containing the input files to use for the input nodes specified
        in pipeline_json.
    remote_project_dir : str
        Path to the project directory on the remote resource.
    output_dir : str
        Path
    id_user : str
        Submitter's id
    pscs_job_id : str
        Job id given by PSCS
    mem_requirement_mb : int
        Amount of memory to request for the job in MB. Default: 5000.
    Returns
    -------
    str
        Path to the generated file
    """
    ht_fd, ht_name = tempfile.mkstemp(prefix="pscs-htcondor-", suffix=".submit")
    # Get HTCondor submit file
    from pscs.static import htcondor_template
    f = open(remapped_input_json, 'r')
    input_files = json.load(f)
    f.close()
    input_paths = ','.join([join(remote_project_dir, str(s)) for s in input_files.values()])

    # need to include other files in input_paths
    input_paths += ',' + join(remote_project_dir, basename(pipeline_json))
    input_paths += ',' + join(remote_project_dir, basename(remapped_input_json))

    start_idx = 0
    if remote_project_dir.startswith(os.sep):
        start_idx = 1
    sep_idx = remote_project_dir.index(os.sep, start_idx)
    output_top = remote_project_dir[:sep_idx]

    # Get log id (new log every 5 days)
    five_days_in_seconds = 5*24*60*60
    current_time = time.time()
    log_label = str(round(current_time - (current_time % five_days_in_seconds)))
    htcondor_proj = htcondor_template.format(node_json=basename(pipeline_json),
                                             input_json=basename(remapped_input_json),
                                             remote_project_dir=remote_project_dir,
                                             input_files=input_paths,
                                             output_directory=output_dir,
                                             output_top=output_top,
                                             log_label=log_label,
                                             pscs_job_id=pscs_job_id,
                                             osf_user_email=current_app.config["REMOTE_COMPUTING_MISC"]["osp"]["NOTIFICATION_EMAIL"],
                                             osf_project_name=current_app.config["REMOTE_COMPUTING_MISC"]["osp"]["PROJECT_NAME"],
                                             sif_path=current_app.config["REMOTE_COMPUTING_MISC"]["osp"]["SIF_PATH"],
                                             mem_requirement_mb=mem_requirement_mb)
    ht_file = os.fdopen(ht_fd, 'w')
    ht_file.write(htcondor_proj)
    ht_file.close()
    return ht_name


def can_user_submit(id_user: str,
                    id_project: str,
                    id_analysis: str,
                    file_info: dict,
                    return_reason: bool = False) -> bool:
    """
    Checks whether a user is allowed to submit another job.
    Parameters
    ----------
    id_user : str
        ID of the user to check.
    id_project : str
        ID of the project for which the job is being submitted.
    id_analysis : str
        Id of the analysis that is being submitted.
    file_info : dict
        Dict of {node_id: {"id" : id_data}} containing the input files to use for the input nodes
    return_reason : bool
        In case the user can't submit, whether a reason should be returned.
    Returns
    -------
    bool
        Whether the user can submit another a job.
    str
        Only returned if return_reason is True. Gives a message explaining why the user can't submit.
    """
    db = get_db()
    # Has the user submitted too many jobs?
    submitted_jobs = get_number_jobs_submitted(id_user, db)
    if submitted_jobs >= current_app.config["SUBMISSION_LIMIT"]:
        if return_reason:
            return False, "Maximum number of jobs already submitted."
        return False

    # Has an identical job already been run?
    # previously_submitted = check_identical_submitted(id_project, id_analysis, file_info, db)
    # if previously_submitted:
    #     if return_reason:
    #         return False, "Identical job has already been submitted."
    #     return False
    if return_reason:
        return True, ""
    else:
        return True


def get_number_jobs_submitted(id_user: str,
                              db: sqlite3.Connection) -> int:
    """
    Returns the number of jobs currently submitted by the user.
    Parameters
    ----------
    id_user : str
        ID of the user to check
    db : sqlite3.Connection
        Database containing job info

    Returns
    -------
    int
        Number of jobs currently submitted
    """
    jobs = db.execute("SELECT id_job "
                      "FROM submitted_jobs "
                      "WHERE id_user = ? AND is_complete = 0", (id_user,)).fetchall()
    return len(jobs)


def check_identical_submitted(id_project: str,
                              id_analysis: str,
                              file_info: dict,
                              db: sqlite3.Connection,
                              ) -> bool:
    """
    Checks whether an identical job has been previously-submitted.
    Parameters
    ----------
    id_project
    id_analysis
    file_info
    db

    Returns
    -------
    bool
        Whether an identical job has been submitted.
    """
    # Convert file_info dict into list of tuples
    file_tuples = []
    for k, v in file_info.items():
        file_tuples.append((k, v))

    previous_jobs = db.execute("SELECT id_job "
                               "FROM submitted_jobs  "
                               "WHERE id_project = ? AND id_analysis = ?", (id_project, id_analysis)).fetchall()
    # check that data + node is all there
    if previous_jobs is None:  # no previous jobs
        return False

    # have previous jobs; check whether data + nodes match
    # every tuple must have a corresponding entry, otherwise the analysis has changed
    for prev_job in previous_jobs:
        for file_tup in file_tuples:
            prev_data = db.execute("SELECT id_job "
                                   "FROM submitted_data "
                                   "WHERE id_job = ? AND id_data = ? AND node_name = ?",
                                   (prev_job["id_job"], file_tup[1], file_tup[0])).fetchall()
            if prev_data is None or len(prev_data) == 0:
                # No match; no need to check the rest
                break
        else:
            # Found match for every data; this is a match
            return True
    return False
