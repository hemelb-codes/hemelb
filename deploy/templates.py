from fabric.api import *
from string import Template
import os

def script_template(template_name):
	source=open(os.path.join(env.localroot,'deploy','templates',template_name))
	result=template(source.read())
	destname=os.path.join(env.localroot,'deploy','.jobscripts',env['name']+'.sh')
	target=open(destname,'w')
	target.write(result)
	return destname
	
def template(pattern):
	return Template(pattern).substitute(env)