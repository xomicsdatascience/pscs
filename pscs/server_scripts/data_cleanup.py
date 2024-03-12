import os
import tempfile
import pathlib
from pathlib import Path
import time
import argparse


tmp_dir = tempfile.gettempdir()
day = 60*60*24
week = day*7
month = day*30
year = day*365.23


def clean_dir(directory: Path,
              age_to_delete_seconds: float):
    """
    Removes files/directories in the specified directory that are older than age_to_delete_seconds.
    Parameters
    ----------
    directory : Path
        Directory in which to perform the cleanup.
    age_to_delete_seconds : float
        Files older than this (in seconds) will be deleted.

    Returns
    -------
    None
    """
    files_to_delete = []
    for file in directory.iterdir():
        if file.is_file() and (time.time() - file.stat().st_mtime) > age_to_delete_seconds:
            files_to_delete.append(file)
    for file in files_to_delete:
        try:
            os.remove(file)
        except PermissionError:
            continue
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", help="Directory in which to perform cleanup", type=str)
    parser.add_argument("--age", help="Age in seconds to delete files", type=float)
    args = parser.parse_args()

    directory = Path(args.directory)
    age_to_delete_seconds = args.age
    clean_dir(directory, age_to_delete_seconds)
