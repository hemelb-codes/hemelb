from templates import *
from fabric.contrib.project import *

@task
def clone():
	run("mkdir -p "+env.remote_directory)
	with cd(env.remote_directory):
		run("rm -rf %s"%env.repository)
		run("hg clone %(hg)s/%(repository)s"%{'hg':env.hg,'repository':env.repository})
	run("mkdir -p %s"%env.build_path)
	run("mkdir -p %s"%env.scripts_path)

@task(alias='cold')
def deploy_cold():
	execute(clone)
	execute(configure)
	execute(build)
	execute(install)
	execute(tools)
	execute(test)

@task
def update_test():
	execute(update)
	execute(configure)
	execute(build)
	execute(install)
	execute(build_python_tools)
	execute(test)

@task
def update():
	with cd(env.repository_path):
		run("hg pull")
		run("hg update")

@task
def clean():
	with cd(env.build_path):
		run("make clean_tree")

@task(alias='tools')
def build_python_tools():
	with cd(env.tools_path):
		run("python setup.py build")

@task
def configure():
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run("cmake %s -DCMAKE_INSTALL_PREFIX=%s -DCMAKE_BUILD_TYPE=%s -DCMAKE_CXX_FLAGS_RELEASE=-O4 -DDEPENDENCIES_INSTALL_PATH=%s" 
			% (env.repository_path, env.install_path, env.build_type, env.install_path)
			)

@task
def build():
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run("make")

@task
def install():
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run("make install")
			run("chmod u+x %s/bin/unittests_hemelb %s/bin/hemelb"%(env.install_path, env.install_path))

@task
def test():
	results_name="test_results.xml"
	with cd(env.remote_directory):
		with prefix(env.run_prefix):
			run(env.pather.join(env.install_path,"bin","unittests_hemelb")+" -o "+results_name)
			get(results_name,os.path.join("remote_files","%(host)s","tests","%(basename)s"))

@task(alias='regress')
def regression_test():
	with cd(env.regression_test_path):
		with prefix(env.python_prefix):
			with prefix(env.run_prefix):
				run("HEMELB_INSTALL_DIR=%s ./diffTest.sh"%env.install_path)
		get("results",os.path.join("remote_files","%(host)s","%(path)s"))
		
@task
def revert(args="--all"):
	with cd(env.repository_path):
		run("hg revert %s"%args)

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
	put("fabric.diff",env.pather.join(env.remote_directory,env.repository))
	with cd(env.repository_path):
		run("patch -p1 < fabric.diff")

@task
def job(template,name,wall_time='0:1:0',nodes='4',memory='1GB'):
	template_name="%s_%s"%(env.template_key,template)
	job_script=fill_in_template(template_name,name=name,
		wall_time=wall_time,nodes=nodes,memory=memory,
		username=env.username,project=env.project,
		executable_path=env.scripts_path
		)
	dest_name=env.pather.join(env.scripts_path,env.pather.basename(job_script))
	put(job_script,dest_name)
	#qsub would go here, but we're just practising for now
	