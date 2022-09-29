"""
This file is intended as a wrapper to abstract away the SQL components of the site. The current development roadmap
has us switching from Python's sqlite3 to a full SQL server later on. To avoid having to go through all the different
calls across the entire codebase, we're centralizing them here.
"""

import sqlite3