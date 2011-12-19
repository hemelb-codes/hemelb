"""
Use python standard library templates to allow strings to include $foo syntax to interpolate elements from the 
Fabric environment dictionary, and to generate job queue submission scripts therefrom.

Job-queue submission scripts should be stored in deploy/templates, with filenames like legion-hemelb (for a script used to
launch hemelb jobs on legion, or hector-unittest, for a script used to launch unit-testing jobs on hector.)
"""
from fabric.api import *
from string import Template
import os

def script_template(template_name):
	"""
	Load a template of the given name, and fill it in based on the Fabric environment dictionary,
	storing the result in deploy/.scripts/job-name.sh
	job-name is loaded from the environment dictionary.
	Return value is the path of the generated script.
	"""
	source=open(os.path.join(env.localroot,'deploy','templates',template_name))
	result=template(source.read())
	destname=os.path.join(env.localroot,'deploy','.jobscripts',env['name']+'.sh')
	target=open(destname,'w')
	target.write(result)
	return destname
	
def template(pattern):
	return Template(pattern).substitute(env)