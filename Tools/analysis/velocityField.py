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
import matplotlib as mpl



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


norm=mpl.colors.Normalize(vmin=0,vmax=0.8)
cmap=cm.jet
m=cm.ScalarMappable(norm=norm, cmap=cmap)



#Iterate through the planes
for child in properties:
    ang=0
    planename = child.attrib.get('file')
    if planename.split('.')[1] == 'xtr' :
        period=int(child.attrib.get('period'))
        #geom=child.find('geometry')
        #normal=geom.find('normal').attrib.get('value')
        #project the normal vector on x,y plane and find a normal vector to it which will be our camera angle
        #normal=normal.replace("(","")
        #normal=normal.replace(")","")
        #normx=float(normal.split(",")[0])
        #normy=float(normal.split(",")[1])
        #viewx=-normy
        #viewy=normx
        #angle=math.acos(viewx/sqrt(viewx**2+viewy**2))     #cos-1(view.x/|view|)
#        if viewy>=0:
#            ax.azim=angle
#        else:
#            ax.azim=-angle

        print planename
        results = VelocityPlaneResults.LoadFromSummary(args.results)
        results.planename=planename


        #Make sure a gif directory exist

        d=res+"/results/Extracted/gifs/vectors_"+planename+"/"
        if not os.path.exists(d):
            os.makedirs(d)

        
        vectors=[]

        #iterate for each period 
        iter_list=[period*(i+1) for i in range(int(lattice/period)-1)]
        for i in iter_list:
        
            if (i*step)<1:    #Ignore time values below 1 sec because our input starts art 1 sec
                continue
            print i
            plane = results.Plane.GetByTimeStep(i)
            xyz = plane.position *pq.metre
            u = plane.velocity * (pq.metre / pq.second)
            tempVel=[]
            x= xyz[:,0]
            y= xyz[:,1]
            z= xyz[:,2]
      
            
            #ax.azim=ang
            #ang=ang+10 
            #if ang >360:
            #    ang=0

            #plot the plane points only for first time
            #if len(vectors)==0:
            p=ax.scatter(x,y,z,c='#e9e9e9',alpha=0.01,s=4)






            modulo=int(len(x)/100)
        

            #plot the velocity vectors
            for index, (px,py,pz,l) in enumerate (zip(x,y,z,u)):
                if index%modulo!=0:
                    continue
                tempVel.append(sqrt(l[0]**2+l[1]**2+l[2]**2))
                ix=float(str(px).split(' ')[0])
                iy=float(str(py).split(' ')[0])
                iz=float(str(pz).split(' ')[0])
                fx=ix+float(str(l[0]).split(' ')[0])/200
                fy=iy+float(str(l[1]).split(' ')[0])/200
                fz=iz+float(str(l[2]).split(' ')[0])/200



                col= m.to_rgba(sqrt(l[0]**2+l[1]**2+l[2]**2))

 
                a = Arrow3D([ix,fx],[iy,fy],[iz,fz], mutation_scale=6, lw=0.7, arrowstyle="-|>", color=col, alpha=0.8)
                vectors.append(ax.add_artist(a))
            plt.axis('off')
            #save the png then delete the scatters and vectors but not the geometry
            plt.savefig(d+str(i)+".png",facecolor='w', edgecolor='w', orientation='portrait',papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)
            #if we have previous vectors remove them
            #if len(vectors)!=0:

            p.remove()
            for vec in vectors:
                vec.remove()
            del vectors[:]     #empty the list of vectors


        file_names = sorted((fn for fn in os.listdir(d) if fn.endswith('.png')))
        images = [Image.open(d+fn) for fn in file_names]
        size=(800,600)
        print "size files="+str(len(file_names))
        print "size images="+str(len(images))
        for im in images:
                
            im.thumbnail(size, Image.ANTIALIAS)
                   
        filename = res+"/results/Extracted/vec_anim_"+planename+".gif"
        writeGif(filename,images, duration=period*step)    


#print epTimestep["position"][i][0], epTimestep["position"][i][1], epTimestep["position"][i][2], epTimestep["velocity"][i][0], epTimestep["velocity"][i][1], epTimestep["velocity"][i][2];

