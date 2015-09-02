
import os
import pyNS
import xml.etree.ElementTree as ET
import os
import re
import math
from os.path import expanduser
import time as ostime
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


#Save the existing results directories in a list to exclude them from the extraction process later
old =[]

for x in os.listdir(results):
    old.append(x)




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
        
        cycles=(4.0/beat)
        cycles=int(cycles)

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
        
        d=devel_lb+"/profiles/Flow_profiles_"+str(i)+"/"
        if not os.path.exists(d):
            os.makedirs(d)


        path = "./Output/Flow/pre/"
        for filename in os.listdir(path):


            filename = os.path.join(path, filename)

            ofilename= filename.split("/")[4]


            
    

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
                for k in range(1,20):     

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
            shutil.copy(ofile,d+ofilename)
        print "Generated velocity output files for "+str(i)


        newsimulation=simulation+"_"+str(i)
        newconfig=devel_lb+"config_files/"+newsimulation+"/"
        if not os.path.exists(newconfig):
            os.makedirs(newconfig)        


        for datafiles in os.listdir(config):
            shutil.copy(config+datafiles,newconfig+datafiles)

        ostime.sleep(10)
#Run simulation using fabric
        os.chdir(deploy)
        os.system("fab " +machine+ " hemelb:config="+newsimulation+",cores="+str(cores)+",wall_time="+wall_time)
        print("\nSumbitted simulation "+str(i)+ "\n\n")
        os.chdir(current_directory)
    f.close()
print("\n\nSIMULATIONS SUBMITTED\n\n")








#Check if all jobs are complete
os.chdir(deploy)
os.system("fab " +machine+ " wait_complete")


#fetch the results
os.system("fab " +machine+ " fetch_results")


print ("\n\nALL RESULTS ARE FETCHED")



#Ignore the old results directories
total=[]
for x in os.listdir(results):
    total.append(x)

s=set(old)

new = [x for x in total if x not in s]


os.chdir(analysis)


for dirname in new:
    print ("Extracting velocities from "+dirname)
    os.system("python allPlanes.py "+results+dirname+"/")
    os.system("python velocityField.py "+results+dirname+"/")
    os.system("python velocityMag.py "+results+dirname+"/")
    os.system("python stress.py "+results+dirname+"/")


os.chdir(current_directory)
print ("ALL EXTRACTIONS DONE")


