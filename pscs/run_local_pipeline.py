import sqlite3
from argparse import ArgumentParser
import os
import subprocess
from pathlib import Path
import os
from uuid import uuid4
from typing import Union
from collections import defaultdict as dd
# from pscs.transfers.fetching import register_results
# from pscs.transfers.dispatching import get_remote_project_dir

# implemented here because there is something happening with imports
# that I can't figure out
image_exts = {"apng", "jpg", "jpeg", "png", "avif", "gif", "svg", "webp", "bmp", "ico", "tiff"}
text_exts = {"txt", "doc", "docx"}
table_exts = {"csv", "tsv"}
data_exts = {"h5ad"}
binary_exts = {"pkl"}

ext_mapper = dd(lambda: "unknown")
for ext in image_exts:
    ext_mapper[ext] = "image"
for ext in text_exts:
    ext_mapper[ext] = "text"
for ext in table_exts:
    ext_mapper[ext] = "table"
for ext in data_exts:
    ext_mapper[ext] = "data"
for ext in binary_exts:
    ext_mapper[ext] = "binary"

def get_file_type(filename: Union[Path, str]) -> str:
    """Examines the extension of the given filename and returns the type of file (e.g., "image", "table", etc.). If
    the type is not recognized, the function returns "unknown"."""
    ext = Path(filename).suffix[1:]
    return ext_mapper[ext]


def register_result(db,
                    id_project,
                    id_analysis,
                    results_directory,
                    interactive_tag = None):
    """Register results into the DB"""
    if isinstance(results_directory, str):
        results_directory = Path(results_directory)
    for result in os.listdir(results_directory):
        file_path = Path(os.path.join(results_directory, result))
        if file_path.name.startswith("pscs-filepaths-"):
            continue
        if not file_path.is_file():
            matching_tag = db.execute("SELECT interactive_tag FROM analysis_interactive_tags WHERE interactive_tag = ?",
                                      (result,)).fetchall()
            if len(matching_tag) == 0:
                continue
            else:
                # Register files in this directory
                register_result(db, id_project, id_analysis, results_directory / result, interactive_tag=result)
                continue
        id_result = str(uuid4())
        file_name = file_path.name
        _, ext = os.path.splitext(file_path)
        new_path = results_directory / (id_result + ext)

        result_type = get_file_type(result)
        # Check if binary; if so, skip; if not, check if binary file exists
        new_binary_path = None
        if result_type == "binary":
            continue
        else:
            binary_file_path = Path(str(file_path) + ".pkl")
            if not binary_file_path.is_file():
                binary_file_path = None
            else:
                new_binary_path = str(new_path) + ".pkl"
                binary_file_path.rename(new_binary_path)
                db.execute("INSERT INTO results_figures "
                           "(id_result, id_result_fig, file_path_fig) VALUES (?,'',?)",
                           (id_result, new_binary_path))
        file_path.rename(new_path)

        # Insert into DB
        if interactive_tag is None:
            db.execute("INSERT INTO results "
                       "(id_result, id_project, id_analysis, file_path, file_name, result_type) "
                       "VALUES (?,?,?,?,?,?) ",
                       (id_result, id_project, id_analysis, str(new_path), file_name, result_type))
        else:
            db.execute("INSERT INTO results "
                       "(id_result, id_project, id_analysis, file_path, file_name, result_type, is_interactive, interactive_tag) "
                       "VALUES (?,?,?,?,?,?,?,?) ",
                       (id_result, id_project, id_analysis, str(new_path), file_name, result_type, 1, interactive_tag))
            if new_binary_path is not None:
                db.execute("INSERT INTO results_figures "
                           "(id_result, id_result_fig, file_path_fig) "
                           "VALUES (?,?,?)",
                           (id_result, id_result, str(new_binary_path)))


    db.commit()
    return


if __name__ == "__main__":
    parser = ArgumentParser(prog = "Local PSCS Pipeline",
                            description="Runs local PSCS pipelines")
    parser.add_argument("db_path", type=str, help="Path to the database file")
    parser.add_argument("instance_url", type=str, help="URL of the Singularity instance")
    parser.add_argument("stdout_file", type=str, help="Path to the stdout file")
    parser.add_argument("stderr_file", type=str, help="Path to the stderr file")
    args = parser.parse_args()
    db_path = args.db_path
    stdout_file = args.stdout_file
    stderr_file = args.stderr_file

    db = sqlite3.connect(db_path)
    db.row_factory = sqlite3.Row
    # Check if there is a currently-running job
    is_currently_running = db.execute("SELECT * FROM local_jobs WHERE is_running = 1").fetchall()
    if len(is_currently_running) > 0:
        exit(0)

    # None running; check job submission
    queued_local = db.execute("SELECT id_job, id_analysis, id_project, remote_results_directory "
                              "FROM submitted_jobs "
                              "WHERE submitted_resource = 'local' AND "
                              "is_complete = 0 "
                              "ORDER BY date_submitted ASC").fetchall()
    while len(queued_local) > 0:
        job_specs = queued_local[0]
        db.execute("UPDATE local_jobs SET is_running = 1 WHERE id_job = ?", (job_specs["id_job"],))
        db.commit()

        # Need: node_json, input_json, output_directory
        pipeline_json = db.execute('SELECT node_file FROM analysis WHERE id_analysis = ?',
                                   (job_specs['id_analysis'],)).fetchone()['node_file']
        proj_dir = os.path.join(job_specs["remote_results_directory"])
        input_json = f"{proj_dir}/pscs-filepaths-{job_specs['id_job']}.json"

        # Run job
        sing_command = f"singularity exec --cleanenv --no-home {args.instance_url} python3.11 /run_pscs_pipeline.py  {pipeline_json} {input_json} {proj_dir}"

        # Mark job as complete
        try:
            subprocess.run(sing_command, shell=True)
            register_result(db,
                            id_project=job_specs["id_project"],
                            id_analysis=job_specs["id_analysis"],
                            results_directory=proj_dir,
                            interactive_tag=None)
        except subprocess.CalledProcessError as e:
            print(f"Error running pipeline: {e}")

        # Record results into DB
        db.execute("UPDATE submitted_jobs "
                   "SET is_complete = 1, "
                   "date_completed=DATETIME('now'), "
                   "stdout_log=?,"
                   "stderr_log=? "
                   "WHERE id_job = ?", (stdout_file, stderr_file, job_specs["id_job"],))
        db.execute("UPDATE local_jobs SET is_running = 0 WHERE id_job = ?", (job_specs["id_job"],))
        db.commit()
        queued_local = db.execute("SELECT id_job, id_analysis, id_project, remote_results_directory "
                                  "FROM submitted_jobs "
                                  "WHERE submitted_resource = 'local' AND "
                                  "is_complete = 0 "
                                  "ORDER BY date_submitted ASC").fetchall()
        if len(queued_local) > 0 and queued_local[0]["id_job"] == job_specs["id_job"]:
            break  # looping