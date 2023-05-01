"""
This file contains code for fetching pipeline results from remote resources. This includes parsing email, fetching data,

"""
from pscs.db import get_db, get_unique_value_for_field
from pscs import parse_env
import subprocess
from pathlib import Path
import os
import sqlite3
from collections import defaultdict as dd
import argparse


def check_for_submitted_jobs(db: sqlite3.Connection):
    """
    Checks whether there are jobs that have been submitted.
    Parameters
    ----------
    db: sqlite3.Connection
        Database to check
    Returns
    -------
    bool
        Whether there are submitted jobs
    """
    incomplete = db.execute("SELECT id_job FROM submitted_jobs WHERE is_complete = 0").fetchone()
    return incomplete is not None


def getmail() -> int:
    """
    Uses the "getmail" utility
    Returns
    -------
    int
        Number of retrieved messages.
    """
    # Note: This is not a stable function; recommend switching to using imapclient.IMAPClient and parsing from there
    response = subprocess.check_output(["getmail", "--new"])
    response_split = response.split("\n")
    num_messages = response_split[-4].split(b' ')[2]
    if num_messages.isnumeric():
        return int(num_messages)
    else:
        return 0


def get_new_mail_files(mail_dir: str, num_new: int) -> list:
    """
    Returns the filenames of the new files
    Parameters
    ----------
    mail_dir : str
        Path to the new mail directory
    num_new : int
        Number of new messages to get.
    Returns
    -------
    list
        List of filepaths pointing to new mail
    """
    mail_list = os.listdir(mail_dir)
    # Order by newest
    mail_list.sort(key=os.path.getctime)
    return mail_list[-num_new:]


def open_mail_file(filepath: str) -> str:
    """Opens the file as an email using the mu utility and returns it"""
    return subprocess.check_output(["mu", "view", filepath])


def get_message_sender(msg: str) -> list:
    """
    Identifies who the message came from
    Parameters
    ----------
    msg : str
        Email message to parse.

    Returns
    -------
    list
        Who sent the message in the form of (user, host)
    """
    # Split by newline
    msg_split = msg.split("\n")
    # Parse until we find "From: "
    for lin in msg_split:
        if lin.startswith("From: "):
            user_host = _extract_between_two_chars(lin, "<", ">")
            return user_host.split("@")
    raise ValueError("No sender found")


def get_job_id(msg: str, resource):
    if resource == "osp":
        return get_osp_job_id(msg)


def get_osp_job_id(msg: str) -> str:
    """
    Extracts the message ID from an email message.
    Parameters
    ----------
    msg : str
        Full email to parse
    Returns
    -------
    str
        OSP job id
    """
    msg_split = msg.split("\n")
    check_string = "Condor job "
    for lin in msg_split:
        if lin.startswith(check_string):
            idx_end = lin.index(".", len(check_string))
            return lin[len(check_string):idx_end]
    return -1


def validate_job_against_db(db: sqlite3.Connection,
                            resource: str,
                            jobid: str) -> bool:
    """
    Checks whether the extracted resource + job id are in the DB and are incomplete
    Parameters
    ----------
    db : sqlite3.Connection
        Database to check
    resource : str
        Resource that was used to obtain results
    jobid : str
        Job ID on the resource

    Returns
    -------
    bool
        Whether the job id is valid and should be fetched.
    """
    # Check db
    is_complete = db.execute("SELECT is_complete "
                             "FROM submitted_jobs "
                             "WHERE submitted_resource = ? AND resource_job_id = ?", (resource, jobid)).fetchone()['is_complete']
    return not is_complete


def fetch_data(remote_info: dict,
               remote_results_dir: str,
               local_results_dir: str,
               ssh_keypath: str):
    """
    Copies the data from the remote computing resource to the local resource.
    Parameters
    ----------
    remote_info : dict
        Info for connecting to the remote server: {"USER": user, "URL": host} as found in the .env file.
    remote_results_dir : str
        Location from which to copy the results.
    local_results_dir
        Location to which to copy the results.
    ssh_keypath : str
        Path to the ssh key to use as identity.
    Returns
    -------
    None
    """
    remote_user = f"{remote_info['USER']}@{remote_info['URL']}"
    out = subprocess.check_output(["rsync", "-re", f"ssh -i{ssh_keypath}", f"{remote_user}:{remote_results_dir}/", local_results_dir])
    return


def update_db_job_complete(db: sqlite3.Connection,
                           resource: str,
                           jobid: str):
    """
    Update the DB to mark the job as completed.
    Parameters
    ----------
    db : sqlite3.Connection
        Database to update.
    resource : str
        Resource on which the results were computed.
    jobid : str
        Remote job id.

    Returns
    -------
    None
    """
    db.execute("UPDATE submitted_jobs "
               "SET is_complete=1, date_completed=DATETIME('now') "
               "WHERE submitted_resource = ? AND resource_job_id = ?", (resource, jobid))
    db.commit()
    return


def _extract_between_two_chars(s: str, char_first: str, char_second: str) -> str:
    """Extracts the contents between of s between char_first and char_second"""
    try:
        idx_first = s.index(char_first)
        idx_second = s.index(char_second, idx_first+1)
        return s[idx_first+1:idx_second]
    except ValueError:
        return ""


def identify_resource(host: str):
    # TODO: use env dictionary to identify resource properly
    if host == "login04.osgconnect.net":
        return "osp"
    else:
        raise ValueError("resource not identifiable")


def poll_resource(resource: str,
                  env: dict) -> dict:
    """
    Checks the specified resource
    Parameters
    ----------
    resource : str
        Resource to poll
    env : dict
        Dict containing the environment info
    Returns
    -------
    dict
        Dict containing job status, keyed by job id
    """
    if resource == "osp":
        return _poll_osp(env)
    else:
        raise NotImplementedError


def _poll_osp(env: dict) -> dict:
    """
    Poll OSP
    Parameters
    ----------
    env : dict
        Dict containing the environment info

    Returns
    -------
    dict
        Dict containing job status, keyed by job id
    """
    # Get condor history
    # TODO: add env to input so we can fetch CRON_KEY
    login_dict = env["REMOTE_COMPUTING"]["osp"]
    login_str = f"{login_dict['USER']}@{login_dict['URL']}"
    cron_key = env["CRON_KEY"]
    history_bin = subprocess.check_output(["ssh", "-o", "IdentitiesOnly=yes", "-i", cron_key, login_str,
                                           "condor_history", "-userlog", "logs/pscs_latest.log"])
    history = history_bin.decode("utf-8")

    history_split = history.split("\n")
    jobs = dict()
    for job_w in history_split[1:]:  # first line is column names
        # remove whitespace
        job = [j for j in job_w.split(' ') if j != ""]
        if len(job) == 0:
            continue
        job_id = job[0].split('.')[0]
        job_status = job[5]
        jobs[job_id] = job_status
    return jobs


def print_debug(s: str, do_print: bool):
    if do_print:
        print(s)
    return


def main(db, env, debug=False):
    # No jobs; exit
    if not check_for_submitted_jobs(db):
        print_debug("No dispatched jobs!", debug)
        return

    # Check messages
    # Poll remote resources
    resources_rows = db.execute("SELECT submitted_resource, resource_job_id "
                                "FROM submitted_jobs "
                                "WHERE is_complete = 0").fetchall()
    resource_set = set()
    resource_job_ids = set()
    num_jobs = dd(int)
    for r in resources_rows:
        res = r["submitted_resource"]
        resource_set.add(res)
        resource_job_ids.add(r["resource_job_id"])
        num_jobs[res] += 1
    print_debug(f"Dispatched jobs: {num_jobs}", debug)
    print_debug(resource_job_ids, debug)
    complete_jobs = []
    print_debug("checking resources", debug)
    for res in resource_set:
        print_debug(f"polling {res}", debug)
        # job_dict = poll_resource(res, env["REMOTE_COMPUTING"][res])
        job_dict = poll_resource(res, env)
        print_debug(job_dict, debug)
        for job_id, job_status in job_dict.items():
            # TODO: "C" works for OSG, not generally
            if job_status == "C":
                complete_jobs.append((res, job_id))
    print_debug(f"completed jobs: {complete_jobs}", debug)
    for c_job in complete_jobs:
        # Transfer data to local dir
        job_info = db.execute("SELECT id_project, id_analysis, remote_results_directory "
                              "FROM submitted_jobs "
                              "WHERE submitted_resource = ? AND resource_job_id = ?", c_job).fetchone()
        if job_info is None:
            print_debug(c_job, debug)
            continue
        print_debug("fetching!", debug)
        results_path = os.path.join(env["INSTANCE_PATH"], "projects", "{id_project}", "results", "{id_analysis}")
        results_path = results_path.format(id_project=job_info["id_project"], id_analysis=job_info["id_analysis"])
        Path(results_path).mkdir(exist_ok=True, parents=True)

        fetch_data(env["REMOTE_COMPUTING"][c_job[0]], job_info["remote_results_directory"], results_path, ssh_keypath=env["CRON_KEY"])
        # register completion
        update_db_job_complete(db, c_job[0], c_job[1])
        for result in os.listdir(results_path):
            print_debug(f"registering {result}", debug)
            id_result = get_unique_value_for_field(db, "id_result", "results")
            file_path = os.path.join("projects", job_info["id_project"], "results", job_info["id_analysis"], result)
            result_type = "result"  # TODO
            db.execute("INSERT INTO results "       
                       "(id_result, id_project, id_analysis, file_path, result_type) "
                       "VALUES (?,?,?,?,?)", (id_result, job_info["id_project"], job_info["id_analysis"], file_path, result_type))
            db.commit()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Results fetcher",
                                     description="Checks the specified database for remotely-submitted jobs; if any, "
                                                 "checks their status on the remote job and copies the results locally "
                                                 "if the jobs are complete.")
    parser.add_argument("--db", help="path to the database file")
    parser.add_argument("--env", help="path to the environment file")
    parser.add_argument("--debug", default=False, action=argparse.BooleanOptionalAction)
    # Daemon is calling; check if there's anything to do
    pargs = vars(parser.parse_args())
    env = parse_env(pargs["env"])
    db = sqlite3.Connection(pargs["db"])
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA foreign_keys = ON;')
    main(db, env, debug=pargs["debug"])

