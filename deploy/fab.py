from templates import *
from machines import *
from fabric.contrib.project import *

@task
def clone():
	run("mkdir -p "+env.remote_path)
	with cd(env.remote_path):
		run("rm -rf %s"%env.repository)
		run("hg clone %(hg)s/%(repository)s"%{'hg':env.hg,'repository':env.repository})
	run("mkdir -p %s"%env.build_path)
	run("mkdir -p %s"%env.scripts_path)
	run("mkdir -p %s"%env.install_path)

@task(alias='cold')
def deploy_cold():
	execute(clone)
	execute(send_distributions)
	execute(configure)
	execute(build)
	execute(install)
	execute(tools)

@task
def update_test():
	execute(update)
	execute(configure)
	execute(build)
	execute(install)
	execute(build_python_tools)

@task
def update():
	with cd(env.repository_path):
		run("hg pull")
		run("hg update")

@task
def clear_build():
	run("rm -rf %s"%env.build_path)
	run("mkdir -p %s"%env.build_path)

@task
def clean():
	with cd(env.build_path):
		run("make clean")

@task(alias='tools')
def build_python_tools():
	with cd(env.tools_path):
		run("python setup.py build")

@task
def configure():
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run("cmake %s -DCMAKE_INSTALL_PREFIX=%s -DDEPENDENCIES_INSTALL_PATH=%s %s" 
			% (env.repository_path, env.install_path, env.install_path, env.cmake_flags)
			)

@task
def build(verbose=False):
	with cd(env.build_path):
		with prefix(env.build_prefix):
			if verbose:
				run("make VERBOSE=1")
			else:
				run("make")

@task
def install():
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run("make install")
			run("chmod u+x %s/bin/unittests_hemelb %s/bin/hemelb"%(env.install_path, env.install_path))

@task
def test():
	with prefix(env.run_prefix):
		execute(job,'unittests','unittests',nodes=1)

@task
def fetch_test_results():
	get(env.pather.join(env.results_path,name,'tests.xml'),
		os.path.join(env.remote_files,"%(host)s","tests","%(basename)s"))

@task(alias='regress')
def regression_test():
	name='regression'
	execute(job,'regression',name)
	get(os.pather.join(env.results_path,name,'results'),
		os.path.join(env.remote_files,"%(host)s","%(path)s"))
		
@task
def revert(args="--all"):
	with cd(env.repository_path):
		run("hg revert %s"%args)
	
@task	
def send_distributions():
	put(os.path.join(env.localroot,'dependencies','distributions','*.tar.gz'),
		env.pather.join(env.repository_path,'dependencies','distributions'))

@task
def sync():
	rsync_project(
		remote_dir=env.repository_path,
		local_dir=env.localroot+'/',
		exclude=map(lambda x: x.replace('\n',''),
		list(open(os.path.join(env.localroot,'.hgignore')))+
		['.hg']+
		list(open(os.path.join(env.localroot,'RegressionTests','.hgignore')))
		)
	)

@task
def patch(args=""):
	local("hg diff %s> fabric.diff"%args)
	put("fabric.diff",env.pather.join(env.remote_path,env.repository))
	with cd(env.repository_path):
		run("patch -p1 < fabric.diff")

@task
def job(template,name,wall_time='0:1:0',nodes=4,memory='1GB'):
	template_name="%s_%s"%(env.machine_name,template)
	results_directory=env.pather.join(env.results_path,name)
	run("mkdir -p %s"%results_directory)
	job_script=fill_in_template(template_name,name=name,
		wall_time=wall_time,nodes=nodes,memory=memory,
		username=env.username,project=env.project,
		executable_path=env.install_path,results=results_directory
		)
	dest_name=env.pather.join(env.scripts_path,env.pather.basename(job_script))
	put(job_script,dest_name)
	run("chmod u+x %s"%dest_name)
	run("%s %s"%(env.job_dispatch,dest_name))
	