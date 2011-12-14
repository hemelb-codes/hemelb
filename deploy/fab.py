from templates import *
from machines import *
from fabric.contrib.project import *
import time

@task
def clone():
	run(template("mkdir -p $remote_path"))
	with cd(env.remote_path):
		run(template("rm -rf $repository"))
		run(template("hg clone $hg/$repository"))

@task(alias='cold')
def deploy_cold():
	execute(clone)
	execute(clear_build)
	execute(prepare_paths)
	execute(send_distributions)
	execute(configure)
	execute(build)
	execute(install)
	execute(build_python_tools)

@task
def update_build():
	execute(update)
	execute(require_recopy)
	execute(configure)
	execute(build)
	execute(install)
	execute(build_python_tools)

@task 
def require_recopy():
	run(template("touch $build_path/hemelb-prefix/src/hemelb-stamp/hemelb-mkdir"))

@task
def update():
	with cd(env.repository_path):
		run("hg pull")
		run("hg update")

@task 
def prepare_paths():
	run(template("mkdir -p $scripts_path"))
	run(template("mkdir -p $results_path"))
	run(template("mkdir -p $config_path"))
	
@task
def clear_build():
	run(template("rm -rf $build_path"))
	run(template("rm -rf $install_path"))
	run(template("mkdir -p $build_path"))
	run(template("mkdir -p $install_path"))
	

@task
def clean():
	with cd(env.build_path):
		run("make clean")

@task(alias='tools')
def build_python_tools():
	with cd(env.tools_path):
		run("python setup.py build")

@task
def stat():
	run("qstat")
	
@task
def monitor():
	while True:
		time.sleep(30)
		execute(stat)

@task
def configure():
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run(template(
			"cmake $repository_path -DCMAKE_INSTALL_PREFIX=$install_path "+
			"-DDEPENDENCIES_INSTALL_PATH=$install_path $cmake_flags"
			))

@task
def build(verbose=False):
	with cd(env.build_path):
		run(template("rm -rf hemelb_prefix/build"))
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
			run(template("chmod u+x $install_path/bin/unittests_hemelb $install_path/bin/hemelb"))
		
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

def with_job(name):
	env.job_results=env.pather.join(env.results_path,name)
	env.job_results_local=os.path.join(env.local_results,name)
	env.job_results_contents=env.pather.join(env.job_results,'*')
	env.job_results_contents_local=os.path.join(env.job_results_local,'*')


def with_config(name):
	env.job_config_path=env.pather.join(env.config_path,name)
	env.job_config_path_local=os.path.join(env.local_configs,name)
	env.job_config_contents=env.pather.join(env.job_config_path,'*')
	env.job_config_contents_local=os.path.join(env.job_config_path_local,'*')
	

@task
def fetch_configs(config=''):
	with_config(config)
	local(template("rsync -pthrvz $username@$remote:$job_config_path/ $job_config_path_local"))

@task
def put_configs(config=''):
	with_config(config)
	run(template("mkdir -p $job_config_path"))
	rsync_project(local_dir=env.job_config_path_local+'/',remote_dir=env.job_config_path)

@task
def put_results(name=''):
	with_job(name)
	run(template("mkdir -p $job_results"))
	rsync_project(local_dir=env.job_results_local+'/',remote_dir=env.job_results)
	
@task
def fetch_results(name=''):
	with_job(name)
	local(template("rsync -pthrvz $username@$remote:$job_results/ $job_results_local"))

@task
def clear_results(name=''):
	with_job(name)
	run(template('rm -rf $job_results_contents'))		

@task
def test():
	execute(job,script='unittests',name='unittests-$build_number-$machine_name',nodes=1)
		
@task
def hemelb(**args):
	options=dict(script='hemelb',
		name='$config-$build_number-$machine_name-$nodes',
		nodes=4,images=10, snapshots=10, steering=1111, wall_time='0:15:0',memory='2G')
	options.update(args)
	execute(put_configs,args['config'])
	execute(job,**options)

@task(alias='regress')
def regression_test():
	execute(job,script='regression',name='regression-$build_number-$machine_name')

@task
def job(**args):
	env.update(name=args['script'],wall_time='0:1:0',nodes=4,memory='1G') #defaults
	env.update(**args)
	env.update(name=template(env.name))
	with_job(env.name)
	
	script_name=template("$machine_name-$script")
	env.job_script=script_template(script_name)
	env.dest_name=env.pather.join(env.scripts_path,env.pather.basename(env.job_script))
	put(env.job_script,env.dest_name)
	
	run(template("mkdir -p $job_results"))
	run(template("cp $dest_name $job_results"))
	run(template("cp $build_cache $job_results"))
	run(template("chmod u+x $dest_name"))
	run(template("$job_dispatch $dest_name"))
	