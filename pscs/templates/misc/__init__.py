"""
This file covers lazy loading of text templates.
"""

import os
from os.path import join
from werkzeug.utils import secure_filename
this_dir = os.path.dirname(__file__)


def __getattr__(name):
    secure_name = secure_filename(name)  # sneaky user
    secure_path = join(this_dir, secure_name)
    secure_path += ".txt"
    if os.path.exists(secure_path):
        f = open(secure_path, 'r')
        dat = f.read()
        f.close()
        return dat
    else:
        raise AttributeError(f"Template {name} was not found.")