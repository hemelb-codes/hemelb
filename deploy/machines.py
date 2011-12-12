from fabric.api import *
import os
import subprocess
import posixpath
import json
from string import Template


env.localroot=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

config=json.load(open(os.path.join(env.localroot,'deploy','machines.json')))
env.update(config)
user_config=json.load(open(os.path.join(env.localroot,'deploy','machines_user.json')))
env.update(user_config)

env.cmake_options={}
env.pather=posixpath
env.job_dispatch=''


@task
def machine(name):
	env.update(config[name])
	env.update(user_config[name])
	env.machine_name=name
	complete_environment()

@task 
def planck():
	execute(machine,'planck')

@task 
def legion():
	execute(machine,'legion')

def complete_environment():
	env.hosts=['%s@%s'%(env.username,env.remote)]
	env.results_path=Template(env.results_path_template).substitute(env)
	env.remote_path=Template(env.remote_path_template).substitute(env)
	env.repository_path=env.pather.join(env.remote_path,env.repository)
	env.tools_path=env.pather.join(env.repository_path,"Tools")
	env.regression_test_path=env.pather.join(env.repository_path,"RegressionTests","diffTest")
	env.tools_build_path=env.pather.join(env.tools_path,'build',env.tools_build)
	env.build_path=env.pather.join(env.remote_path,'build')
	env.install_path=env.pather.join(env.remote_path,'install')
	env.scripts_path=env.pather.join(env.remote_path,'scripts')
	env.local_results=os.path.join(env.localroot,'results')
	env.cmake_total_options=env.cmake_default_options.copy()
	env.cmake_total_options.update(env.cmake_options)
	env.cmake_flags=' '.join(["-D%s=%s"%option for option in env.cmake_total_options.iteritems()])
	module_commands=["module %s"%module for module in env.modules]
	env.build_prefix=" && ".join(module_commands+env.build_prefix_commands) or 'echo Building...'
	env.run_prefix=" && ".join(module_commands+env.run_prefix_commands) or 'echo Running...'
	env.python_prefix="export PYTHONPATH=$PYTHONPATH:%s"%env.pather.join(env.tools_build_path)
	env.build_number=subprocess.check_output(['hg','id','-i','-rtip','%s/%s'%(env.hg,env.repository)]).strip()
	