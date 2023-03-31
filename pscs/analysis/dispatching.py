# This file contains code for dispatching processing pipelines to different computational resources.
import sqlite3
import subprocess
from pscs.db import get_db
from flask import current_app
from os.path import join, basename
import json
import tempfile
import os


def dispatch(pipeline_json: str,
             file_ids: dict,
             id_project: str,
             id_analysis: str,
             resource: str = 'osp'):
    """
    Collects the relevant files and sends them to the relevant computational resource.
    Parameters
    ----------
    pipeline_json : str
        Path to the JSON file containing the node information.
    file_ids : dict
        Dict of {node_id: file_path} containing the input files to use for the input nodes specified in pipeline_json.
    id_project : str
        ID of the project for which the analysis is being run.
    id_analysis : str
        ID of the analysis being run.
    resource : str
        Name of the resource that should be used.
    Returns
    -------
    None
    """
    # Make destination
    remote_proj_dir = get_remote_project_dir(id_project=id_project, id_analysis=id_analysis, resource=resource)
    remote_mkdir(remote_dir=remote_proj_dir, resource=resource)
    remote_results = join(remote_proj_dir, 'results')
    remote_mkdir(remote_dir=remote_results, resource=resource)

    # Transfer files
    transfer_file(pipeline_json, remote_dir=remote_proj_dir, resource=resource)
    transfer_files(local_files=list(file_ids.values()), remote_dir=remote_proj_dir, resource=resource)

    # Create and transfer dictionary for new filepaths on remote machine
    remote_paths = remap_filepaths(file_ids, remote_dir=remote_proj_dir)
    remap_filename = save_filepath_remap(remote_paths)
    transfer_file(remap_filename, remote_dir=remote_proj_dir, resource=resource)

    # Submission script
    htcondor_script = generate_htcondor_submission(pipeline_json=pipeline_json,
                                                   remapped_input_json=remap_filename,
                                                   remote_project_dir=remote_proj_dir,
                                                   output_dir=remote_results)
    transfer_file(htcondor_script, remote_proj_dir, resource=resource)
    cmd = f"condor_submit {os.path.join(remote_proj_dir, basename(htcondor_script))}"
    remote_cmd(cmd)
    return


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


def remap_filepaths(file_ids: dict,
                    remote_dir: str) -> dict:
    """
    Converts the input file paths to paths on the remote resource.
    Parameters
    ----------
    file_ids : dict
        Dict of {node_id: file_path} containing the input files to use for the input nodes specified in pipeline_json.
    remote_dir : str
        Path on the remote resource.

    Returns
    -------
    dict
        Dict of {node_id: file_path} with the file_path updated to reflect the path on the remote resource.
    """
    remote_paths = {}
    for node_id, file_path in file_ids.items():
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
                                 ) -> str:
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
    input_paths = ','.join([join(remote_project_dir,str(s)) for s in input_files.values()])

    # need to include other files in input_paths
    input_paths += ',' + join(remote_project_dir, basename(pipeline_json))
    input_paths += ',' + join(remote_project_dir, basename(remapped_input_json))
    input_paths += ',' + "test.py"

    start_idx = 0
    if remote_project_dir.startswith(os.sep):
        start_idx = 1
    sep_idx = remote_project_dir.index(os.sep, start_idx)
    output_top = remote_project_dir[:sep_idx]

    htcondor_proj = htcondor_template.format(node_json=basename(pipeline_json),
                                             input_json=basename(remapped_input_json),
                                             remote_project_dir=remote_project_dir,
                                             input_files=input_paths,
                                             output_directory=output_dir,
                                             output_top=output_top)
    ht_file = os.fdopen(ht_fd, 'w')
    ht_file.write(htcondor_proj)
    ht_file.close()
    return ht_name
