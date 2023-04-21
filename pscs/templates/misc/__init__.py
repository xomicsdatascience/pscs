import os
from os.path import join
this_dir = os.path.dirname(__file__)

f = open(join(this_dir, 'confirmation_email.txt'))
confirmation_template = f.read()
f.close()