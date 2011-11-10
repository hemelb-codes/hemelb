from machines import *

@task
def clone():
	run("mkdir -p "+env.remote_directory)
	with cd(env.remote_directory):
		run("rm -rf %s"%env.repository)
		run("hg clone %(hg)s/%(repository)s"%{'hg':env.hg,'repository':env.repository})
	
@task(alias='cold')
def deploy_cold():
	execute(clone)
	execute(parmetis)
	execute(build)
	execute(build_python_tools)
	execute(test)

@task
def update():
	with cd(env.repository_path):
		run("hg pull")
		run("hg update")

@task
def clean():
	with cd(env.code_path):
		run("make clean_tree")

@task
def configure_parmetis():
	with cd(env.parmetis_path):
		with prefix(env.build_prefix):
			run("make config")

@task
def build_parmetis():
	with cd(env.parmetis_path):
		with prefix(env.build_prefix):
				run("make")

@task(alias="metis")
def parmetis():
	execute(configure_parmetis)
	execute(build_parmetis)

@task(alias='tools')
def build_python_tools():
	with cd(env.tools_path):
		run("python setup.py build")

@task
def build():
	with cd(env.code_path):
		with prefix(env.build_prefix):
			run("make")
		
@task
def test():
	results_name="test_results.xml"
	with cd(env.remote_directory):
		run(env.pather.join(env.code_path,"build","unitTests")+" 2>"+results_name)
		get(results_name,os.path.join("remote_files","%(host)s","tests","%(basename)s"))

@task(alias='regress')
def regression_test():
	with cd(env.regression_test_path):
		with prefix(env.python_prefix):
			run("./diffTest.sh")
		get("results",os.path.join("remote_files","%(host)s","%(path)s"))
		
@task
def revert(args="--all"):
	with cd(env.repository_path):
		run("hg revert %s"%args)
		
@task
def patch(args=""):
	local("hg diff %s> fabric.diff"%args)
	put("fabric.diff",env.pather.join(env.remote_directory,env.repository))
	with cd(env.repository_path):
		run("patch -p1 < fabric.diff")