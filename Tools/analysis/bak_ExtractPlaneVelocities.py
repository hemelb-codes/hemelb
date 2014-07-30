#!/usr/bin/env python
import os.path
from hemeTools.parsers.extraction import ExtractedProperty
from analysis.HemeLbRunResults import HemeLbRunResults
import quantities as pq
from math import sqrt
#from matplotlib import cm
#from pylab import plot,show,contour
#from scipy.cluster.vq import whiten
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib
#import matplotlib.pyplot as plt
import pdb
import xml.etree.ElementTree as ET


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

    #parse the arguments which are (the result directory)
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('results', help='results directory to process')
    
    args = p.parse_args()
    

    
    res=args.results
    res=str(res)

    #get the step size, the plane names and period from the config file
    tree = ET.parse(res+'/config.xml')
    root = tree.getroot()
    step = root[0][0].attrib.get('value')
    step= float(step)
    lattice = root[0][1].attrib.get('value')
    lattice = int(lattice)

    #Iterate through the defined planes in xml and extraxt the corresponding xtr file
    properties = root.find('properties')
    for child in properties:
        planename = child.attrib.get('file')
        if planename.split('.')[1] == 'xtr' and planename!='planeo1.xtr':     #planeo1 was excluded because it is a corrupted file
            fname=res+"/results/Extracted/"+planename.split('.')[0]+".txt"
            fh1=open(fname,'w')
            results = VelocityPlaneResults.LoadFromSummary(args.results)
            results.planename=planename
            period=int(child.attrib.get('period'))
            
            #for each period without reaching the total number of lattice
            iter_list=[period*(i+1) for i in range(int(lattice/period))]  #int(lattice/period)  in place of the 3
            for i in iter_list:
                print i
                
                plane = results.Plane.GetByTimeStep(i)
                xyz = plane.position *pq.metre
                u = plane.velocity * (pq.metre / pq.second)
                tempVel=[]
                # Compute the maximum velocity in the plane perpendicular to the centreline
                for l in u:
                    tempVel.append(sqrt(l[0]**2+l[1]**2+l[2]**2))
                k=map(str,[i*step])[0]
                out_txt=k+' '+map(str,[max(tempVel)])[0]
                print >> fh1, out_txt
                
            
        
    

