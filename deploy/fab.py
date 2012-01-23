"""
Fabric definitions for HemeLB
Usage:
        With current working directory anywhere inside the HemeLB mercurial checkout, execute
                "fab <machinename> <task>"
        for example:
                "fab hector deploy_cold"

        Do fab -l to get a list of available commands.

        Before use, you MUST copy deploy/machines_user_example.yml as deploy/machines_user.yml and fill in your personal details.
        For smoothest usage, you should also ensure ssh keys are in place so the target machine can access the mercurial repository
        of hemelb.
"""

from templates import *
from machines import *
from fabric.contrib.project import *
import time
import re
import numpy as np

@task
def clone():
    """Delete and checkout the repository afresh."""
    run(template("mkdir -p $remote_path"))
    if env.no_ssh or env.no_hg:
        with cd(env.remote_path):
            run(template("rm -rf $repository"))
        # Some machines do not allow outgoing connections back to the mercurial server
        # so the data must be sent by a project sync instead.
        execute(sync)
         # On such machines, we cannot rely on an outgoing connection to servers to find dependencies either.
    else:
        with cd(env.remote_path):
            run(template("rm -rf $repository"))
            run(template("hg clone $hg/$repository"))
    if env.no_ssh or env.needs_tarballs:
        execute(send_distributions)
    execute(copy_regression_tests)

@task(alias='cold')
def deploy_cold():
    """Checkout, build, and install hemelb, from new."""
    execute(clear_build)
    execute(clone)
    execute(prepare_paths)
    execute(configure)
    execute(build)
    execute(install)
    execute(build_python_tools)

@task
def update_build():
    """Update and do incremental build."""
    execute(update)
    execute(require_recopy)
    execute(configure)
    execute(build)
    execute(install)
    execute(build_python_tools)

@task
def require_recopy():
    """Notify the build system that the code has changed."""
    run(template("touch $build_path/hemelb-prefix/src/hemelb-stamp/hemelb-mkdir"))
    run(template("rm -rf $build_path/hemelb-prefix"))

@task
def update():
    """Update the remote mercurial repository"""
    if env.no_ssh:
        execute(sync)
    else:
        with cd(env.repository_path):
            run("hg pull")
            run("hg update")

@task
def prepare_paths():
    """Create remote locations to store results and configs"""
    run(template("mkdir -p $scripts_path"))
    run(template("mkdir -p $results_path"))
    run(template("mkdir -p $config_path"))
    run(template("mkdir -p $profiles_path"))

@task
def clear_build():
    """Wipe out existing built and installed HemeLB."""
    run(template("rm -rf $build_path"))
    run(template("rm -rf $install_path"))
    run(template("mkdir -p $build_path"))
    run(template("mkdir -p $install_path"))
    run(template("mkdir -p $temp_path"))


@task
def clean():
    """Clean remote build area."""
    with cd(env.build_path):
        run("make clean")

@task(alias='tools')
def build_python_tools():
    """Build and install python scripts."""
    with cd(env.tools_path):
        with prefix(env.build_prefix):
            run(template("python setup.py install --prefix $install_path"))

@task
def stat():
    """Check the remote message queue status"""
    #TODO: Respect varying remote machine queue systems.
    run(template("qstat -u $username"))

@task
def monitor():
    """Report on the queue status, ctrl-C to interrupt"""
    while True:
        execute(stat)
        time.sleep(30)

@task
def configure():
    """CMake configure step for HemeLB and dependencies."""
    with cd(env.build_path):
        with prefix(env.build_prefix):
            run(template("rm -f $build_path/CMakeCache.txt"))
            run(template(
            "cmake $repository_path -DCMAKE_INSTALL_PREFIX=$install_path "+
            "-DHEMELB_DEPENDENCIES_INSTALL_PATH=$install_path $cmake_flags "+
            "-DHEMELB_SUBPROJECT_MAKE_JOBS=$make_jobs"
            ))

@task
def code_only():
    """Configure, build, and install for the /Code C++ code only, do not attempt to install and build dependencies."""
    execute(configure_code_only)
    execute(build_code_only)
    execute(install_code_only)

@task
def configure_code_only():
    """CMake configure step for HemeLB code only."""
    run(template("rm -rf $code_build_path"))
    run(template("mkdir -p $code_build_path"))
    with cd(env.code_build_path):
        with prefix(env.build_prefix):
            run(template(
            "cmake $repository_path/Code -DCMAKE_INSTALL_PREFIX=$install_path "+
            "-DHEMELB_DEPENDENCIES_INSTALL_PATH=$install_path $cmake_flags"
            ))

@task
def build_code_only(verbose=False):
    """CMake build step for HemeLB code only."""
    with cd(env.code_build_path):
        with prefix(env.build_prefix):
            if verbose or env.verbose:
                run(template("make -j$make_jobs VERBOSE=1"))
            else:
                run(template("make -j$make_jobs"))

@task
def install_code_only():
    """CMake install step for HemeLB code only."""
    with cd(env.code_build_path):
        with prefix(env.build_prefix):
            run("make install")
            run(template("chmod u+x $install_path/bin/unittests_hemelb $install_path/bin/hemelb"))

@task
def build(verbose=False):
    """CMake build step for HemeLB and dependencies."""
    with cd(env.build_path):
        run(template("rm -rf hemelb_prefix/build"))
        with prefix(env.build_prefix):
            if verbose or env.verbose:
                run("make VERBOSE=1")
            else:
                run("make")

@task
def install():
    """CMake install step for HemeLB and dependencies."""
    with cd(env.build_path):
        with prefix(env.build_prefix):
            #run("make install") // Doesn't have a separate install step, as HemeLB gets installed by the sub-project to the final install dir.
            run(template("chmod u+x $install_path/bin/unittests_hemelb $install_path/bin/hemelb"))

@task
def revert(args="--all"):
    """Revert local changes in the remote repository.
     Including those made through 'fab ... patch' and 'fab ... sync'.
    Specify a path relative to the repository root to revert only some files or directories
    """
    with cd(env.repository_path):
        run("hg revert %s"%args)

@task
def send_distributions():
    """Transmit dependency tarballs to remote.
    Files taken from dependencies/distributions.
    Useful to prepare a build on target machines with CMake before 2.8.4, where
    HTTP redirects are not followed.
    """
    run(template("mkdir -p $repository_path/dependencies/distributions"))
    rsync_project(local_dir=os.path.join(env.localroot,'dependencies','distributions')+'/',
    remote_dir=env.pather.join(env.repository_path,'dependencies','distributions'))

@task
def fetch_distributions():
    """Fetch dependency tarballs tfrom remote.
    Files taken from dependencies/distributions.
    Useful to prepare a build on target machines with CMake before 2.8.4, where
    HTTP redirects are not followed.
    """
    local(template("rsync -pthrvz $username@$remote:$repository_path/dependencies/distributions/ $localroot/dependencies/distributions"))

@task
def copy_regression_tests():
    if env.regression_test_source_path != env.regression_test_path:
        run(template("cp -r $regression_test_source_path $regression_test_path"))

@task
def sync():
    """Update the remote repository with local changes.
    Uses rysnc.
    Respects the local .hgignore files to avoid sending unnecessary information.
    """
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
    """Update the remote repository with local changes.
    Uses hg diff to generate patchfiles.
    Specify a path relative to the repository root to patch only some files or directories
    e.g 'fab legion patch Code/main.cc'
    """
    local("hg diff %s> fabric.diff"%args)
    put("fabric.diff",env.pather.join(env.remote_path,env.repository))
    with cd(env.repository_path):
        run("patch -p1 < fabric.diff")

def with_job(name):
    """Augment the fabric environment with information regarding a particular job name.
    Definitions created:
    job_results: the remote location where job results should be stored
    job_results_local: the local location where job results should be stored
    """
    env.job_results=env.pather.join(env.results_path,name)
    env.job_results_local=os.path.join(env.local_results,name)
    env.job_results_contents=env.pather.join(env.job_results,'*')
    env.job_results_contents_local=os.path.join(env.job_results_local,'*')


def with_config(name):
    """Internal: augment the fabric environment with information regarding a particular configuration name.
    Definitions created:
    job_config_path: the remote location where the config files for the job should be stored
    job_config_path_local: the local location where the config files for the job may be found
    """
    env.job_config_path=env.pather.join(env.config_path,name)
    env.job_config_path_local=os.path.join(env.local_configs,name)
    env.job_config_contents=env.pather.join(env.job_config_path,'*')
    env.job_config_contents_local=os.path.join(env.job_config_path_local,'*')

def with_profile(name):
    """Internal: augment the fabric environment with information regarding a particular profile name.
    Definitions created:
    job_profile_path: the remote location where the profile should be stored
    job_profile_path_local: the local location where the profile files may be found
    """
    env.profile=name
    env.job_profile_path=env.pather.join(env.profiles_path,name)
    env.job_profile_path_local=os.path.join(env.local_profiles,name)
    env.job_profile_contents=env.pather.join(env.job_profile_path,'*')
    env.job_profile_contents_local=os.path.join(env.job_profile_path_local,'*')

@task
def fetch_configs(config=''):
    """
    Fetch config files from the remote, via rsync.
    Specify a config directory, such as 'cylinder' to copy just one config.
    Config files are stored as, e.g. cylinder/config.dat and cylinder/config.xml
    Local path to use is specified in machines_user.json, and should normally point to a mount on entropy,
    i.e. /store4/blood/username/config_files
    This method is not intended for normal use, but is useful when the local machine cannot have an entropy mount,
    so that files can be copied to a local machine from entropy, and then transferred to the compute machine,
    via 'fab entropy fetch_configs; fab legion put_configs'
    """
    with_config(config)
    local(template("rsync -pthrvz $username@$remote:$job_config_path/ $job_config_path_local"))

@task
def put_configs(config=''):
    """
    Transfer config files to the remote.
    For use in launching jobs, via rsync.
    Specify a config directory, such as 'cylinder' to copy just one configuration.
    Config files are stored as, e.g. cylinder/config.dat and cylinder/config.xml
    Local path to find config directories is specified in machines_user.json, and should normally point to a mount on entropy,
    i.e. /store4/blood/username/config_files
    If you can't mount entropy, 'fetch_configs' can be useful, via 'fab entropy fetch_configs; fab legion put_configs'
    """
    with_config(config)
    run(template("mkdir -p $job_config_path"))
    rsync_project(local_dir=env.job_config_path_local+'/',remote_dir=env.job_config_path)

@task
def put_results(name=''):
    """
    Transfer result files to a remote.
    Local path to find result directories is specified in machines_user.json.
    This method is not intended for normal use, but is useful when the local machine cannot have an entropy mount,
    so that results from a local machine can be sent to entropy, via 'fab legion fetch_results; fab entropy put_results'
    """
    with_job(name)
    run(template("mkdir -p $job_results"))
    rsync_project(local_dir=env.job_results_local+'/',remote_dir=env.job_results)

@task
def fetch_results(name=''):
    """
    Fetch results of remote jobs to local result store.
    Specify a job name to transfer just one job.
    Local path to store results is specified in machines_user.json, and should normally point to a mount on entropy,
    i.e. /store4/blood/username/results.
    If you can't mount entropy, 'put results' can be useful,  via 'fab legion fetch_results; fab entropy put_results'
    """
    with_job(name)
    local(template("rsync -pthrvz $username@$remote:$job_results/ $job_results_local"))

@task
def clear_results(name=''):
    """Completely wipe all result files from the remote."""
    with_job(name)
    run(template('rm -rf $job_results_contents'))

@task
def fetch_profiles(name=''):
    """
    Fetch results of remote jobs to local result store.
    Specify a job name to transfer just one job.
    Local path to store results is specified in machines_user.json, and should normally point to a mount on entropy,
    i.e. /store4/blood/username/results.
    If you can't mount entropy, 'put results' can be useful,  via 'fab legion fetch_results; fab entropy put_results'
    """
    with_profile(name)
    local(template("rsync -pthrvz $username@$remote:$job_profile_path/ $job_profile_path_local"))

@task
def put_profiles(name=''):
    """
    Transfer result files to a remote.
    Local path to find result directories is specified in machines_user.json.
    This method is not intended for normal use, but is useful when the local machine cannot have an entropy mount,
    so that results from a local machine can be sent to entropy, via 'fab legion fetch_results; fab entropy put_results'
    """
    with_profile(name)
    run(template("mkdir -p $job_profile_path"))
    rsync_project(local_dir=env.job_profile_path_local+'/',remote_dir=env.job_profile_path)

@task(alias='test')
def unit_test(**args):
    """Submit a unit-testing job to the remote queue."""
    options=dict(script='unittests',name='unittests_${build_number}_${machine_name}',cores=1)
    options.update(args)
    execute(job,**options)

@task
def hemelb(**args):
    """Submit a HemeLB job to the remote queue.
    The job results will be stored with a name pattern of $config_${build_number}_${machine_name}_$cores,
    e.g. cylinder-abcd1234-legion-256
    Keyword arguments:
            config : config directory to use to define geometry, e.g. config=cylinder
            cores : number of compute cores to request
            images : number of images to take
            snapshots : number of snapshots to take
            steering : steering session i.d.
            wall_time : wall-time job limit
            memory : memory per node
    """
    options=dict(script='hemelb',
            name='$config_${build_number}_${machine_name}_$cores',
            cores=4,images=10, snapshots=10, steering=1111, wall_time='0:15:0',memory='2G')
    options.update(args)
    execute(put_configs,args['config'])
    execute(job,**options)

@task(alias='regress')
def regression_test(**args):
    """Submit a regression-testing job to the remote queue."""
    execute(copy_regression_tests)
    options=dict(name='regression_${build_number}_${machine_name}',cores=3,
            wall_time='0:20:0', images=1, snapshots=1, steering=1111,script='regression')
    options.update(args)
    execute(job,**options)

@task
def job(**args):
    """Internal low level job launcher.
    Execute a generic job on the remote machine. Use hemelb, regress, or test instead."""
    env.update(name=args['script'],wall_time='0:1:0',cores=4,memory='1G') #defaults
    env.update(**args)
    env.update(name=template(env.name))
    # If we request less than one node's worth of cores, need to keep N<=n
    if int(env.corespernode)>int(env.cores):
        env.corespernode=env.cores
    env['job_name']=env.name[0:env.max_job_name_chars]
    with_job(env.name)

    script_name=template("$template_key-$script")
    env.job_script=script_template(script_name)

    env.dest_name=env.pather.join(env.scripts_path,env.pather.basename(env.job_script))
    put(env.job_script,env.dest_name)

    run(template("mkdir -p $job_results"))
    run(template("cp $dest_name $job_results"))
    run(template("cp $build_cache $job_results"))
    run(template("chmod u+x $dest_name"))
    with cd(env.job_results):
        run(template("$job_dispatch $dest_name"))

def input_to_range(arg,default):
    ttype=type(default)
    gen_regexp="\[([\d\.]+):([\d\.]+):([\d\.]+)\]" #regexp for a array generator like [1.2:3:0.2]
    if not arg:
        return [default]
    match=re.match(gen_regexp,arg)
    if match:
        vals=list(map(ttype,match.groups()))
        if ttype==int:
            return range(*vals)
        else:
            return np.arange(*vals)

    return [ttype(arg)]

@task
def create_config(profile,config_template="${profile}_${VoxelSize}_${Steps}_${Cycles}",VoxelSize=None,Steps=None,Cycles=None):
    """Create a config file
    Create a config file (geometry and xml) based on a configuration profile.
    """
    with_profile(profile)

    from HemeLbSetupTool.Model.Profile import Profile
    p = Profile()
    p.LoadFromFile(os.path.expanduser(os.path.join(env.job_profile_path_local,profile)+'.pro'))
    p.StlFile=os.path.expanduser(os.path.join(env.job_profile_path_local,profile)+'.stl')
    env.VoxelSize=str(VoxelSize).replace(".","_")
    env.Steps=Steps
    env.Cycles=Cycles
    config=template(config_template)
    with_config(config)
    local(template("mkdir -p $job_config_path_local"))
    p.OutputConfigFile=os.path.expanduser(os.path.join(env.job_config_path_local,'config.dat'))
    p.OutputXmlFile=os.path.expanduser(os.path.join(env.job_config_path_local,'config.xml'))
    p.VoxelSize=VoxelSize
    p.Steps=Steps
    p.Cycles=Cycles
    if not p.IsReadyToGenerate:
        raise "Not ready to generate"
    p.Generate()

@task
def create_configs(profile,config_template="${profile}_${VoxelSize}_${Steps}_${Cycles}",VoxelSize=None,Steps=None,Cycles=None):
    """Create many config files, by looping over multiple voxel sizes, step counts, or cycles
    """
    with_profile(profile)
    from HemeLbSetupTool.Model.Profile import Profile
    p = Profile()
    p.LoadFromFile(os.path.expanduser(os.path.join(env.job_profile_path_local,profile)+'.pro'))
    for currentVoxelSize in input_to_range(VoxelSize,p.VoxelSize):
        for currentSteps in input_to_range(Steps,1000):
            for currentCycles in input_to_range(Cycles,3):
                execute(create_config,profile,config_template,currentVoxelSize,currentSteps,currentCycles)
                
@task 
def hemelb_profile(profile,config_template="${profile}_${VoxelSize}_${Steps}_${Cycles}",VoxelSize=None,Steps=None,Cycles=None,**args):
    """Submit HemeLB job(s) to the remote queue.
    The HemeLB config file(s) will be prepared according to the profile and profile arguments given.
    This can submit a massive number of jobs -- do not use on systems with a limit to number of queued jobs permitted.
    """
    with_profile(profile)
    from HemeLbSetupTool.Model.Profile import Profile
    p = Profile()
    p.LoadFromFile(os.path.expanduser(os.path.join(env.job_profile_path_local,profile)+'.pro'))
    for currentVoxelSize in input_to_range(VoxelSize,p.VoxelSize):
        for currentSteps in input_to_range(Steps,1000):
            for currentCycles in input_to_range(Cycles,3):
                for currentCores in input_to_range(args['cores'],4):
                    execute(create_config,profile,config_template,currentVoxelSize,currentSteps,currentCycles)
                    hemeconfig={}
                    hemeconfig.update(args)
                    hemeconfig['config']=template(config_template)
                    hemeconfig['cores']=currentCores
                    execute(hemelb,hemeconfig)