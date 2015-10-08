
import os
import pyNS
import xml.etree.ElementTree as ET
import os
import re
import math
from os.path import expanduser
import time
import shutil


#PAth to home
home = expanduser("~")



#You may need to change the path to the directories below

#Path to results directory
results = home+"/devel/ccs/fabric/hemelb/results/"

#Path to HemeLB Tools/analysis
analysis= home+"/hemelb-dev/hemelb/Tools/analysis/"


#Path to pyNS-master
current_directory=os.path.dirname(os.path.realpath(__file__))



#Path to hemelb deploy
deploy = home+"/hemelb-dev/hemelb/deploy"








#parse the arguments which are (the number of cores)
import argparse
p = argparse.ArgumentParser()
p.add_argument('machine', help='machine name')
p.add_argument('cores', help='number of cores')
p.add_argument('simulation', help='name of config')
p.add_argument('wall_time', help='wall_time for each simulation')
args = p.parse_args()


machine=args.machine


#Number of cores
cores = int(args.cores)
simulation=args.simulation
wall_time=args.wall_time



#Path to the config file
config = home+"/devel/ccs/fabric/hemelb/config_files/"+simulation+"/"



#Path to the hemelb in devel
devel_lb = home+"/devel/ccs/fabric/hemelb/"





profiles=devel_lb+"profiles/"

i=0

for direc in os.listdir(profiles):
    i=i+1

    newsimulation=simulation+"_"+str(i)
    newconfig=devel_lb+"config_files/"+newsimulation+"/"
    if not os.path.exists(newconfig):
        os.makedirs(newconfig)        

    for filename in os.listdir(profiles+direc):
        shutil.copy(profiles+direc+"/"+filename,config)

    for datafiles in os.listdir(config):
        shutil.copy(config+datafiles,newconfig+datafiles)

    ostime.sleep(10)


    
    os.chdir(deploy)
    os.system("fab " +machine+ " hemelb:config="+newsimulation+",cores="+str(cores)+",wall_time="+wall_time)
    os.chdir(current_directory)
    for filename in os.listdir(profiles+direc):
        os.remove(config+filename)
    time.sleep(5)
    print "Submitted job "+str(i)


print "\n\nALL JOBS ARE SUBMITTED\n\n"



