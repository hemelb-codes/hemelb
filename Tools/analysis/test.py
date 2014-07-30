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



if __name__ == '__main__':


    

    hemeMagic = 0x686C6221 # hlb!
    extractMagic = 0x78747204 # xtr
    geometryMagic = 0x676D7904 # gmy number


    



    velocity=True
    stress=True
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
   



           
    #Make directories for png files
    d=res+"/results/Extracted/gifs/display/"
    if not os.path.exists(d):
        os.makedirs(d)


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




    period=30000

    properties = root.find('properties')

    figs=[]

    print "\n\nExtracting plane velocities\n\n"
    #Iterate through the defined planes in xml and extraxt the corresponding xtr file
    iter_list=[period*(i+1) for i in range(int(lattice/period)-1)]
    for i in iter_list:
        if (i*step)<1:    #Ignore time values below 1 sec because our input starts art 1 sec
            continue
        print i    

        
    
        for child in properties:
            planename = child.attrib.get('file')
            if planename=="plane49.xtr" or planename=="plane59.xtr":     
                results = VelocityPlaneResults.LoadFromSummary(args.results)
                results.planename=planename
                
                plane = results.Plane.GetByTimeStep(i)
                xyz = plane.position *pq.metre
                u = plane.velocity * (pq.metre / pq.second)
                tempVel=[]
                for l in u:
                    tempVel.append(sqrt(l[0]**2+l[1]**2+l[2]**2))
                #Plot the points velocities of the plane

                x= xyz[:,0]
                y= xyz[:,1]
                z= xyz[:,2]
                
                figs.append(ax.scatter(x,y,z,c=tempVel,cmap=cm.jet, s=3))

                figs[len(figs)-1].set_clim(vmin=0,vmax=1)
                

        print len(figs)
        #plt.colorbar(figs[0])                
        ax.azim=-14
        plt.axis('off')
        #plt.title(str(i*step)+ " sec")

        plt.savefig(d+str(i)+".png",facecolor='w', edgecolor='w', orientation='portrait',papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)
        
        for fig in figs:
            fig.remove()
        del figs[:]     #empty the list of figures   

#Convert the png files to gif animation
file_names = sorted((fn for fn in os.listdir(d) if fn.endswith('.png')))
images = [Image.open(d+fn) for fn in file_names]
size=(800,600)
print "size files="+str(len(file_names))
print "size images="+str(len(images))
for im in images:
                
    im.thumbnail(size, Image.ANTIALIAS)
                #im.convert('RGB').convert('P',palette=Image.ADAPTIVE)
                
filename = res+"/results/Extracted/display.gif"
writeGif(filename,images, duration=period*step)    




            


