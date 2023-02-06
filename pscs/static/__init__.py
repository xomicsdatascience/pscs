import pkg_resources

file = pkg_resources.resource_filename('pscs', 'static/submission_htcondor.submit')
f = open(file, 'r')
htcondor_template = f.read()
f.close()
