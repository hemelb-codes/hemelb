### RSDT-hosted Jenkins server


RSDT hosts a jenkins server at http://development.rc.ucl.ac.uk/jenkins/. The hemelb jobs defined
there are available to the heme-team for perusal. Login is done using *github* credentials. Please
email rc-softdev@ucl.ac.uk with *github* credentials to add new users.


The jobs are generated automatically from the a set of yaml files living in a [dedicated orphan
branch on github](https://github.com/UCL-CCS/hemelb-dev/tree/test_jobs). The grammar for the yaml
files is given by the [jenkins-job-builder](http://ci.openstack.org/jenkins-job-builder/). There is
a specific job, called setup-jenkins-jobs-heme that can be built to (re-)generate the jobs described
in the folder.  *Alternatively*, the jobs can be updated from a Heme repository as follows:


    # install python-jenkinsm, jenkins-job-builder and two plugins for the latter
    > pip install python-jenkins==0.4.8 # newer version seem to bug out on Mac
    > pip install jenkins-job-builder
    > pip install git+https://github.com/UCL/jenkjobs.git
    > pip install jenkins-jobs-slack
    # add ucl-ccs-slack-token from slack integration
    # found here, originally: https://ucl-ccs.slack.com/services/B0BT7QA10
    > echo "token from slack integration" > /path/to/heme/source/ucl-ccs-slack-token
    # create/update jobs
    > jenkins-jobs --conf /path/to/configfile.ini update /path/to/heme-repo/ [some-job-heme]

The configuration file tells jenkins-job-builder where to find and how to authenticate with the
jenkins server. It goes something like this:

     [jenkins]
     user=mdavezac
     password= some *jenkins* token
     url=http://development.rc.ucl.ac.uk/jenkins/


The token can be obtained from jenkins under people/you/configure/show api token. Or email
rc-softdev@ucl.ac.uk for help. The jenkins token (to access jenkins) and the slack token (for
jenkins to ping slack) are different API tojens. Both are necessary.

- to add dedicated tests for long-lived branches, edit hemelb.yml and copy/paste a new project
- to modify the build, look at the builder in template.yml.

Some caveats:

- all jobs must end with heme or hemelb
- jobs are under version control, so using the jenkins gui results in temporary changes only
