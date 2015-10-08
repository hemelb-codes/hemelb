#!/usr/bin/env python
import os.path
from hemeTools.parsers.extraction import ExtractedProperty
from analysis.HemeLbRunResults import HemeLbRunResults
import quantities as pq
from math import sqrt
from matplotlib import cm
#from pylab import plot,show,contour
#from scipy.cluster.vq import whiten
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib
import pdb
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import matplotlib.animation as animation
from PIL import Image
from images2gif import writeGif
import math



hemeMagic = 0x686C6221 # hlb!
extractMagic = 0x78747204 # xtr
geometryMagic = 0x676D7904 # gmy



class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)







class VelocityPlaneResults(HemeLbRunResults):
    @classmethod
    def LoadFromSummary(cls, run):
        inst = cls()
        inst.RunDirectory = run
        return inst


def ComputeDistanceToLine(points, linePoint, lineVector):
    """ Computes the distance between points and a line defined by the linePoint
        point and the lineVector vector
    """
    vectorFromLineToPoints = (points - linePoint) - np.dot((points - linePoint), lineVector)[:, np.newaxis] * lineVector
    return np.apply_along_axis(np.linalg.norm, -1, vectorFromLineToPoints) * pq.metre



#parse the arguments which are (the result directory)
import argparse
p = argparse.ArgumentParser()
p.add_argument('results', help='results directory to process')
    
args = p.parse_args()
    
res=args.results
res=str(res)



#get the step size, the plane names and period from the config file
tree = ET.parse(res+'config.xml')
root = tree.getroot()
step = root[0][0].attrib.get('value')
step= float(step)
lattice = root[0][1].attrib.get('value')
lattice = int(lattice)

datperiod=0
datname=""

properties = root.find('properties')
for child in properties:
    filename = child.attrib.get('file')
    if filename.split('.')[1] == 'dat':
        datperiod=int(child.attrib.get('period'))
        datname=filename



epSet = ExtractedProperty(res+"/results/Extracted/"+datname);
epTimestep = epSet.GetByTimeStep(2*datperiod);

#The geometry plotting


x=[]
y=[]
z=[]

for i in range(1,epSet.siteCount):
    if str(epTimestep["d_shearstress"][i])!="inf":
        print i
        x.append(epTimestep["position"][i][0])
        y.append(epTimestep["position"][i][1])
        z.append(epTimestep["position"][i][2])
 
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

g=ax.scatter(x,y,z,c='#e9e9e9',s=4,alpha=0.01)






#Iterate through the planes
for child in properties:
    ang=0
    planename = child.attrib.get('file')
    if planename.split('.')[1] == 'xtr' :

        period=int(child.attrib.get('period'))

        results = VelocityPlaneResults.LoadFromSummary(args.results)
        results.planename=planename

        print planename, period
        plane = results.Plane.GetByTimeStep(period)
        xyz = plane.position *pq.metre

        x= xyz[:,0]
        y= xyz[:,1]
        z= xyz[:,2]
        index=np.argmax(z)
        ax.scatter(x,y,z,c='b',alpha=0.1,s=4)
        
        maxx=float(str(x[index]).split(' ')[0])
        maxy=float(str(y[index]).split(' ')[0])
        maxz=float(str(z[index]).split(' ')[0])
        maxx,maxy, _ = proj3d.proj_transform(maxx,maxy,maxz,ax.get_proj())
        ax.annotate(planename.split('.')[0], xycoords='data',xy=(maxx,maxy),xytext=(0,50),textcoords='offset points', ha= 'right', va='bottom',bbox=dict(boxstyle='round,pad=0.5',fc='yellow',alpha=0.5),arrowprops= dict(arrowstyle= '->',connectionstyle='arc3,rad=0'))


plt.axis('off')

        
            #save the png then delete the scatters and vectors but not the geometry
plt.savefig(res+"/results/Extracted/allPlanes.png",facecolor='w', edgecolor='w', orientation='portrait',papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)
            #if we have previous vectors remove them
            #if len(vectors)!=0:

