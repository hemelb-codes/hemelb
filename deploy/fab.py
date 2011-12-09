from machines import *

@task
def clone():
	run("mkdir -p "+env.remote_directory)
	with cd(env.remote_directory):
		run("rm -rf %s"%env.repository)
		run("hg clone %(hg)s/%(repository)s"%{'hg':env.hg,'repository':env.repository})
		run("mkdir -p dependencies/build")
		run("mkdir -p Code/build")
	
@task(alias='cold')
def deploy_cold():
	execute(clone)
	execute(configure)
	execute(build)
	execute(install)
	execute(tools)
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



@task(alias='tools')
def build_python_tools():
	with cd(env.tools_path):
		run("python setup.py build")

@task
def configure():
	with cd(env.code_build_path):
		with prefix(env.build_prefix):
			run("ccmake .. -DCMAKE_INSTALL_PREFIX=%s"%env.install_path)

@task
def build():
	with cd(env.code_build_path):
		with prefix(env.build_prefix):
			run("make")

def install():
	with cd(env.code_build_path):
		with prefix(env.build_prefix):
			run("make install")

@task
def test():
	results_name="test_results.xml"
	with cd(env.remote_directory):
		run(env.pather.join(env.install_path,"bin","unitTests")+" 2>"+results_name)
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