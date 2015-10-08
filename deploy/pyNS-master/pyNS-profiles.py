
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
p.add_argument('simulation', help='name of config')

args = p.parse_args()



#Number of cores

simulation=args.simulation




#Path to the config file
config = home+"/devel/ccs/fabric/hemelb/config_files/"+simulation+"/"



#Path to the hemelb in devel
devel_lb = home+"/devel/ccs/fabric/hemelb/"







#Get values from InputSet.csv
import csv
with open ('InputSet.csv','rb') as f:
    reader =csv.reader(f)
    data=list(reader)
    lines=len(data)
        
    for i in range(1,lines):
            
        meanPressure=data[i][0]
        cardiacOutput=data[i][1]
        dynamicViscosity=data[i][2]
        heartRate=data[i][3]

        hrate=float(heartRate)
        beat=60.0/hrate
        Period=str(beat)
#Modify the boundary_conditions_willis.xml with the values extracted above
        tree = ET.parse('XML/Models/WillisCircle/boundary_conditions_willis.xml')
        root = tree.getroot()

        root[1][0][0][0].text=meanPressure
        root[1][0][1][0].text=cardiacOutput
        root[1][0][3][0].text=dynamicViscosity
        root[1][0][2].text=Period

        tree.write('XML/Models/WillisCircle/boundary_conditions_willis.xml')

        print("\nRunning pyNS for data "+str(i)+"\n\n")

#Run pyNS for the Circle of Willis
        try:
            pyNS.runSimulation("specific", False, 'XML/', 'Output/', 'Images/', 'XML/XSD/', None, 'vascular_mesh_v1.1.xml', 'vascular_output.xml', None, 'vascular_network_v3.2.xsd', 'boundary_conditions_v3.1.xsd', 'willis', 'XML/parameters.csv', False, -1, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, None)
        except KeyboardInterrupt:
            print "\nLocal web server for post processing was shutdown successfully. pyNS is ready for next simulation."
            print "If you want to inspect last simulation results, type ./pyNS.py --results last"
            print "If you want to inspect a specific simulation results, type ./pyNS.py --results folderName"

#######################################################


#Take the flow output files and generate velocities output files
#Each element in pyNS has two radiuses: Minimum and Maximum
#Assumption: We will take the average radius to get the output velocities
        
        d=devel_lb+"pyNS/Flow_profiles_"+str(i)+"/"
        if not os.path.exists(d):
            os.makedirs(d)


        path = "./Output/Flow/pre/"
        for filename in os.listdir(path):


            filename = os.path.join(path, filename)

            ofilename= filename.split("/")[4]


            shutil.copy(filename,d+ofilename)
    

        print "Generated flow profiles for "+str(i)

    f.close()
print("\n\nALL FLOW PROFILES ARE GENERATED\n\n")



