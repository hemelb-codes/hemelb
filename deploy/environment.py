from fabric.api import *
import os
import posixpath
env.hg='ssh://hg@entropy.chem.ucl.ac.uk'
env.repository='hemelb'
env.modules=[]
env.pather=posixpath
env.remote_directory='FabricHemeLb'

def complete_environment():
	env.repository_path=env.pather.join(env.remote_directory,env.repository)
	env.tools_path=env.pather.join(env.repository_path,"Tools")
	env.regression_test_path=env.pather.join(env.repository_path,"RegressionTests","diffTest")
	env.tools_build_path=env.pather.join(env.tools_path,'build',env.tools_build)
	env.build_path=env.pather.join(env.remote_directory,'build')
	env.install_path=env.pather.join(env.remote_directory,'install')
	module_commands=["module %s"%module for module in env.modules]
	env.build_prefix=" && ".join(module_commands+["export HEMELB_MACHINE=%s"%env.HEMELB_MACHINE])
	env.python_prefix="export PYTHONPATH=$PYTHONPATH:%s"%env.pather.join("~",env.tools_build_path)