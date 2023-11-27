from flask import Blueprint, current_app
import os
import shutil
bp = Blueprint('homepage', __name__, url_prefix="/")


def get_posts(paths: list,
              check_cache: bool = True,
              write_cache: bool = True):
    """
    Loads data from the files listed in paths. If check_cache is True, checks whether the files are stored in the
    tmpfs directory. If write_cache is True and the file isn't found in tmpfs, also copies over the file to tmpfs.
    Parameters
    ----------
    paths : list
        List of paths to files containing the post information.
    check_cache : bool
        Whether to check the tmpfs cache first.
    write_cache : bool
        Whether to copy the file to tmpfs if it isn't already there.

    Returns
    -------
    list
        List of dicts containing the post information.
    """
    cache_dir = os.path.join(current_app.config["TMPFS"], "posts")
    posts = []
    for p in paths:
        filename = os.path.filename(p)
        if check_cache:
            cache_file = os.path.join(cache_dir, filename)
            if os.path.exists(cache_file):
                # Data exists and is cached
                posts.append(read_post(cache_file))
                continue
        if os.path.exists(p):
            # Data exists
            posts.append(read_post(p))
            if write_cache:
                copy_to_cache(p, "posts")
    return posts


def read_post(path: str):
    """Reads the text file specified by path and returns it."""
    f = open(path, "r")
    contents = f.read()
    f.close()
    return contents


def copy_to_cache(filepath: str,
                  cache_dest: str = ""):
    """
    Copies the specified path to tmpfs.
    Parameters
    ----------
    filepath : str
        Path of the file to cache.
    cache_dest : str
        Name of the directory under which to store files.

    Returns
    -------
    None
    """
    cache_dir = os.path.join(current_app.config["TMPFS"], cache_dest)
    filename = os.path.filename(filepath)
    shutil.copy2(filepath, os.path.join(cache_dir, filename))
    return
