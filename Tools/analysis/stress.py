#!/usr/bin/env python
import os
import xdrlib
import numpy as np
from hemeTools.parsers.extraction import ExtractedProperty
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from PIL import Image
from images2gif import writeGif
from mpl_toolkits.mplot3d import Axes3D
from math import sqrt
from matplotlib import cm
import xml.etree.ElementTree as ET


hemeMagic = 0x686C6221 # hlb!
extractMagic = 0x78747204 # xtr
geometryMagic = 0x676D7904 # gmy



#parse the arguments which are (the result directory)
import argparse
p = argparse.ArgumentParser()
p.add_argument('results', help='results directory to process')
    
args = p.parse_args()
    
res=args.results
res=str(res)




#Make directories for png files
d=res+"/results/Extracted/gifs/stress/"
if not os.path.exists(d):
    os.makedirs(d)





#get the step size, the plane names and period from the config file
tree = ET.parse(res+'config.xml')
root = tree.getroot()
step = root[0][0].attrib.get('value')
step= float(step)
lattice = root[0][1].attrib.get('value')
lattice = int(lattice)


print "\n\nExtracting maximum wall shear stress\n\n"
#Get the name of the datafile (.dat)
properties = root.find('properties')
for child in properties:
    filename = child.attrib.get('file')
    if filename.split('.')[1] == 'dat':
        datafile=filename
        period=int(child.attrib.get('period'))

epSet = ExtractedProperty(res+"/results/Extracted/"+datafile);
maxstress=[]

time=[]
ang=0

write_ascii = True

iter_list=[period*(i+1) for i in range(int(lattice/period)-1)]
for i in iter_list:
    if (i*step)<1:    #Ignore time values below 1 sec because our input starts art 1 sec
        continue

    x=[]
    y=[]
    z=[]    


    text_file = None
    if(write_ascii):
      text_file = open(d+datafile+"."+str(i)+".ascii.txt", "w")

    print i
    epTimestep = epSet.GetByTimeStep(i);

    stress=[]
    for j in range(1,epSet.siteCount):
        if str(epTimestep["d_shearstress"][j])!="inf":
            stress.append(epTimestep["d_shearstress"][j])
            x.append(epTimestep["position"][j][0])
            y.append(epTimestep["position"][j][1])
            z.append(epTimestep["position"][j][2])

            if(write_ascii):
              text_file.write(str(epTimestep["position"][j][0])+" "+str(epTimestep["position"][j][1])+" "+str(epTimestep["position"][j][2])+" "+str(epTimestep["d_shearstress"][j])+"\n")
    
    maxstress.append(max(stress))
    print "max stress = "+str(max(stress))
    print "min stress = "+str(min(stress))
    time.append(i*step)
    
    #TEMP: exit after writing ASCII files!
    if(write_ascii):
        continue

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    g=ax.scatter(x,y,z,c=stress,cmap=cm.jet,s=4,alpha=0.4)
    ax.azim=ang
    plt.colorbar(g)
    g.set_clim(vmin=0,vmax=20)
    ang=ang-1
    plt.title(str(i*step)+ " sec")
    plt.axis('off')
    plt.savefig(d+str(i)+".png",facecolor='w', edgecolor='w', orientation='portrait',papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

    plt.cla()    
    plt.clf()


#Convert the png files to gif animation
file_names = sorted((fn for fn in os.listdir(d) if fn.endswith('.png')))
images = [Image.open(d+fn) for fn in file_names]
size=(800,600)
print "size files="+str(len(file_names))
print "size images="+str(len(images))
for im in images:
               
    im.thumbnail(size, Image.ANTIALIAS)
                #im.convert('RGB').convert('P',palette=Image.ADAPTIVE)
                
filename = res+"/results/Extracted/Stress.gif"
writeGif(filename,images, duration=period*step*2)    
   


#Plot max shear stress
plt.plot(time,maxstress)
plt.axis('on')
plt.xlabel('time (s)')
plt.ylabel('Max Stress (N/m2)')
            
plt.savefig(res+"/results/Extracted/MaxStress.png",facecolor='w', edgecolor='w', orientation='portrait',papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)
plt.clf()
plt.cla()   



