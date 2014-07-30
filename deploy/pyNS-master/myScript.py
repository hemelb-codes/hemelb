
import os
import pyNS
import xml.etree.ElementTree as ET
import os
import re
import math
from os.path import expanduser


#Path to pyNS-master
current_directory=os.path.dirname(os.path.realpath(__file__))

#PAth to home
home = expanduser("~")

#Path to hemelb deploy
deploy = home+"/hemelb-dev/hemelb/deploy"

#Simulation name (name of config directory)
simulation = "nhnn_mar2014"

#Path to the config file
config = home+"/devel/ccs/fabric/hemelb/config_files/"+simulation+"/"

#Number of cores
cores = 3






#Get values from InputData.csv
import csv
with open ('InputSet.csv','rb') as f:
        reader =csv.reader(f)
        data=list(reader)
        lines=len(data)
        
        for i in range(1,lines):
            
            meanPressure=data[i][0]
            cardiacOutput=data[i][1]
            dynamicViscosity=data[i][2]

#Modify the boundary_conditions_willis.xml with the values extracted above
            tree = ET.parse('XML/Models/WillisCircle/boundary_conditions_willis.xml')
            root = tree.getroot()

            root[1][0][0][0].text=meanPressure
            root[1][0][1][0].text=cardiacOutput
            root[1][0][3][0].text=dynamicViscosity

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

##            dirname = "Velocities"
##            if not os.path.isdir("./" + dirname + "/"):
##                os.mkdir("./" + dirname + "/")
##
##            subdirname = "Velocities-"+str(i)
##            if not os.path.isdir("./Velocities/" + subdirname + "/"):
##                os.mkdir("./Velocities/" + subdirname + "/")
##

            path = "./Output/Flow/pre/"
            for filename in os.listdir(path):
                filename = os.path.join(path, filename)
                ofilename= filename.split("/")[4]
    
#                ofile= os.path.join("./Velocities/" + subdirname + "/",ofilename)
                ofile = os.path.join(config,ofilename)
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
                    for line in file_object:
                        flow=float(line)
                        vel=(flow/(math.pi*rad*rad))*((10**-6)/60)
            
                        velofile.write(str(time)+"\t"+str(vel) + "\n")
                        time=time+0.005
                    velofile.write("2.0\t0")
                    velofile.close()

            print "Generated velocity output files for "+str(i)


#Run simulation using fabric
            os.chdir(deploy)
            os.system("fab oppenheimer hemelb:config="+simulation+",cores="+str(cores))
            print("\nSumbitted simulation "+str(i)+ "\n\n")
            os.chdir(current_directory)
        f.close()
print("\n\nSIMULATIONS SUBMITTED\n\n")
print("Use stat and fetch_results to get the results")
