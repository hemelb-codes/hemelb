from machines import *
from string import Template

def fill_in_template(template_name,**binding):
	source=open(os.path.join(env.localroot,'deploy','templates',template_name))
	result=Template(source.read()).substitute(binding)
	destname=os.path.join(env.localroot,'deploy','templates',"%s.sh"%binding['name'])
	target=open(destname,'w')
	target.write(result)
	return destname