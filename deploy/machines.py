from environment import *
	
@task
def planck():
	env.template_key='legion' #not really, but for testing
	env.username='jamespjh'
	env.project='jamespjh_ccs'
	env.hosts=['planck.chem.ucl.ac.uk']
	env.tools_build="lib.linux-x86_64-2.6"
	complete_environment()

@task
def legion():
	env.template_key='legion'
	env.username='ucgajhe'
	env.project='jamespjh_ccs'
	env.hosts=["%s@legion.rc.ucl.ac.uk"%env.username]
	env.tools_build="lib.linux-x86_64-2.6"
	env.modules.append("remove ofed/qlogic/intel/64/1.2.7")
	env.modules.append("add ofed/openmpi/gcc/64/1.4.1")
	env.modules.append("load cmake")
	complete_environment()