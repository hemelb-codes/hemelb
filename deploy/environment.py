from fabric.api import *
import os
import posixpath
env.hg='ssh://hg@entropy.chem.ucl.ac.uk'
env.repository='hemelb'
env.modules=[]
env.pather=posixpath
env.remote_directory_name='FabricHemeLb'
env.build_type='Release'

def complete_environment():
	env.remote_directory=env.pather.join("~",env.remote_directory_name)
	env.repository_path=env.pather.join(env.remote_directory,env.repository)
	env.tools_path=env.pather.join(env.repository_path,"Tools")
	env.regression_test_path=env.pather.join(env.repository_path,"RegressionTests","diffTest")
	env.tools_build_path=env.pather.join(env.tools_path,'build',env.tools_build)
	env.build_path=env.pather.join(env.remote_directory,'build')
	env.install_path=env.pather.join(env.remote_directory,'install')
	module_commands=["module %s"%module for module in env.modules]
	env.build_prefix=" && ".join(module_commands)
	env.run_prefix=env.build_prefix + ("&& LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s/dependencies/lib"%env.repository_path)
	env.python_prefix="export PYTHONPATH=$PYTHONPATH:%s"%env.pather.join(env.tools_build_path)