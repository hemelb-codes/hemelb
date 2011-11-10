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
	#env.modules.append("purge")
	#env.modules.append("load ucl")
	env.modules.append("load cmake")
	#env.modules.append("load gcc")
	env.modules.append("unload ofed/qlogic")
	env.modules.append("load ofed/openmpi/intel/64/1.4.1")
	complete_environment()