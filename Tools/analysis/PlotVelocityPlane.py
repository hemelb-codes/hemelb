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
# Example of use: env PYTHONPATH=../../irregular-boundaries/scripts/:$PYTHONPATH python PlotVelocityPlane.py ~/FabricHemeLb/results/newsurface_aneuFine_hector_080543757d99+_20131104123855/

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
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('results', help='results directory to process')
    args = p.parse_args()
    fname=os.path.dirname(os.path.abspath(__file__))+'/velocity_profile_59_C.txt'
    fh1=open(fname,'w')
    results = VelocityPlaneResults.LoadFromSummary(args.results)
    iter_list=[100000*(i+1) for i in range(1)]

    fname = 'plane00.xtr'
    results.planename=fname
    for i in iter_list:
        print i
        
        plane = results.Plane.GetByTimeStep(i)
        xyz = plane.position *pq.metre
        u = plane.velocity * (pq.metre / pq.second)
        tempVel=[]
        # Compute the maximum velocity in the plane perpendicular to the centreline
        for l in u:
            tempVel.append(sqrt(l[0]**2+l[1]**2+l[2]**2))
        k=map(str,[i*3.02e-06])[0]
        out_txt=k+' '+map(str,[max(tempVel)])[0]
        print >> fh1, out_txt
    # Use matplotlib to do a 3D plot of the velocity magnitude in a plane
#    fig = plt.figure()
#    x = plane.position[:,0]
#    y = plane.position[:,1]
#    z = plane.position[:,2]
#    X,Y=np.meshgrid(x,y)
#    pdb.set_trace()
#    ax = fig.add_subplot(111,projection='3d')
#    tv=np.array(tempVel)
#    fe = plane.position
    #fe=np.loadtxt(open("/home/jitta/position.csv","rb"),dtype='float64',delimiter=",")
#    f=whiten(fe)
#    x1 = f[:,0]
#    y1 = f[:,1]
#    z1 = f[:,2]
#    ax.scatter3D(x1,z1,tv)
#    show()
