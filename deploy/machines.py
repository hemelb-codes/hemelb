from fabric.api import *
import os
import posixpath
env.update(
	hg='ssh://hg@entropy.chem.ucl.ac.uk',
	repository='hemelb',
	modules=[],
	pather=posixpath,
	remote_directory_name='FabricHemeLb',
	build_type='Release',
	localroot=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
	build_prefix_commands=[],
	run_prefix_commands=[]
)

machines={
	'legion':dict(
		username='ucgajhe',
		project='jamespjh_ccs',
		job_dispatch='qsub',
		remote="legion.rc.ucl.ac.uk",
		tools_build="lib.linux-x86_64-2.6",
		modules=["remove ofed/qlogic/intel/64/1.2.7", "add ofed/openmpi/gcc/64/1.4.1", "load cmake"],
		results_path='/home/ucgajhe/Scratch/FabricHemeLb/results',
		remote_path='/home/ucgajhe/FabricHemeLb'
	),
	'planck':dict(
		remote='planck.chem.ucl.ac.uk',
		job_dispatch='',
		username='jamespjh',
		project='jamespjh_ccs',
		hosts=['planck.chem.ucl.ac.uk'],
		tools_build="lib.linux-x86_64-2.6",
		results_path='/home/jamespjh/FabricHemeLb/results',
		remote_path='/home/jamespjh/FabricHemeLb'
	)
}

@task
def machine(name):
	env.update(machines[name])
	env.machine_name=name
	complete_environment()

def complete_environment():
	env.hosts=['%s@%s'%(env.username,env.remote)]
	env.repository_path=env.pather.join(env.remote_path,env.repository)
	env.tools_path=env.pather.join(env.repository_path,"Tools")
	env.regression_test_path=env.pather.join(env.repository_path,"RegressionTests","diffTest")
	env.tools_build_path=env.pather.join(env.tools_path,'build',env.tools_build)
	env.build_path=env.pather.join(env.remote_path,'build')
	env.install_path=env.pather.join(env.remote_path,'install')
	env.scripts_path=env.pather.join(env.remote_path,'scripts')
	env.remote_files=os.path.join(env.localroot,'remote_files')
	module_commands=["module %s"%module for module in env.modules]
	env.build_prefix=" && ".join(module_commands+env.build_prefix_commands) or 'echo Building...'
	env.run_prefix=" && ".join(module_commands+env.run_prefix_commands) or 'echo Running...'
	env.python_prefix="export PYTHONPATH=$PYTHONPATH:%s"%env.pather.join(env.tools_build_path)
	