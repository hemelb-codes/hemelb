from environment import *
	
@task
def planck():
	env.hosts=['planck']
	env.HEMELB_MACHINE='CCS'
	env.tools_build="lib.linux-x86_64-2.6"
	complete_environment()

@task
def legion():
	env.hosts=['ucgajhe@legion.rc.ucl.ac.uk']
	env.HEMELB_MACHINE='Legion'
	env.tools_build="lib.linux-x86_64-2.6"
	env.modules.append("remove ofed/qlogic/intel/64/1.2.7")
	env.modules.append("add ofed/openmpi/gcc/64/1.4.1")
	env.modules.append("load cmake")
	complete_environment()