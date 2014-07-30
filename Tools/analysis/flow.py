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

d=res+"/results/Extracted/gifs/"
if not os.path.exists(d):
    os.makedirs(d)


epSet = ExtractedProperty(res+"/results/Extracted/flow_snapshot.dat");
#print epSet.voxelSizeMetres
#print epSet.originMetres
#print epSet.siteCount
#print epSet.fieldCount
#print epSet.times
#print epSet.GetFieldSpec()


iter_list=[300*(i+1) for i in range(int(22000/300))]  #0.75 was added since the simulation did not run for the whole time steps
for i in iter_list:
    if (i*300)<1:    #Ignore time values below 1 sec because our input starts art 1 sec
        continue
    print i
    epTimestep = epSet.GetByTimeStep(i);
    x=[]
    y=[]
    z=[] 
    vel=[]               

    for j in range(1,epSet.siteCount):
        
        x.append(epTimestep["position"][j][0])
        y.append(epTimestep["position"][j][1])
        z.append(epTimestep["position"][j][2])
        vel.append(sqrt(epTimestep["d_velocity"][j][0]**2+epTimestep["d_velocity"][j][1]**2+epTimestep["d_velocity"][j][2]**2))


                             

    #Plot the points velocities of the plane
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

                
    g=ax.scatter(x,y,z,c=vel,cmap=cm.OrRd)


    plt.colorbar(g)
    g.set_clim(vmin=0,vmax=0.1)
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
                
filename = res+"/results/Extracted/velo_animation.gif"
writeGif(filename,images, duration=300*9.09*(10**(-5)))




#for i in range(1,epSet.siteCount):
#    print epTimestep["position"][i][0], epTimestep["position"][i][1], epTimestep["position"][i][2],epTimestep["shearstress"][i];

#for i in range(1,epSet.siteCount):
#    print epTimestep["position"][i][0], epTimestep["position"][i][1], epTimestep["position"][i][2],epTimestep["pressure"][i];

#print epTimestep["position"][i][0], epTimestep["position"][i][1], epTimestep["position"][i][2], epTimestep["velocity"][i][0], epTimestep["velocity"][i][1], epTimestep["velocity"][i][2];

