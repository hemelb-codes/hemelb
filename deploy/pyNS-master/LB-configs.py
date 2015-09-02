
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


#Path to the hemelb in devel
devel_lb = home+"/devel/ccs/fabric/hemelb/"







#parse the arguments which are (the number of cores)
import argparse
p = argparse.ArgumentParser()
p.add_argument('simulation', help='name of config')

args = p.parse_args()



#Number of cores

simulation=args.simulation




#Path to the config file
config = home+"/devel/ccs/fabric/hemelb/config_files/"+simulation+"/"


i=0

flowprofiles=devel_lb+"pyNS/"

for profile in os.listdir(flowprofiles):
    i=i+1
    d=devel_lb+"profiles/Velocity_profiles_"+str(i)+"/"
    if not os.path.exists(d):
        os.makedirs(d)

    flowdir=flowprofiles+profile+"/"

    for filename in os.listdir(flowdir):
        ofilename= filename
        filename = os.path.join(flowdir, filename)

        ofile = os.path.join(d,ofilename)
        velofile = open(ofile, 'w')
        velofile.write("0.0\t0\n")
        with open(filename) as file_object:
            line = file_object.readline()
        
            line = file_object.readline()
        
            minrad=line.split(" ")[2]
            line = file_object.readline()
            maxrad=line.split(" ")[2]
            minrad = float(minrad)*0.001
            maxrad = float(maxrad)*0.001
            rad=(minrad+maxrad)/2
            time=1.0
            for k in range(1,10):      #3 cardiac cycle  between 1 sec and 4 secs

                for line in file_object:
                    flow=float(line)
                    vel=(flow/(math.pi*rad*rad))*((10**-6)/60)

                    velofile.write(str(time)+"\t"+str(vel) + "\n")
                    time=time+0.005
                    if time>3.995:
                        break

                if time>3.995:
                    break
                
                file_object.seek(0)
                line = file_object.readline()
                line = file_object.readline()
                line = file_object.readline()
            velofile.write("4.0\t0")
            velofile.close()
      
    print "Generated velocity config files for "+str(i)



print "\n\nALL LB config velocity files are generated\n\n"




