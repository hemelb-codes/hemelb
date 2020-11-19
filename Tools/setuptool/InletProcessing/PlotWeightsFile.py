# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys


#420.500000 297.500000 215.000000

def filterLowX(xyzv):
    deleteList = []
    for i in xrange(0, len(xyzv)):
        if xyzv[i,3] > 0.99:
            deleteList.append(i)
            continue
        if xyzv[i,0] > 423 or xyzv[i,0] < 418:
            deleteList.append(i)
            continue
        if xyzv[i,1] > 300 or xyzv[i,1] < 295:
            deleteList.append(i)
            continue
        if xyzv[i,2] > 219 or xyzv[i,2] < 212:
            deleteList.append(i)
            continue

    deleteList.sort(reverse=True)

    for i in deleteList:
        xyzv = np.delete(xyzv,i,0)

    return xyzv



if __name__ == "__main__":

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('v')

    for i in xrange(1,len(sys.argv)):
        xyzv = np.loadtxt(sys.argv[i]) 
        #xyzv = filterLowX(np.loadtxt(sys.argv[i]))
      
        mode = "xyv"
        cl='blue'
        if(i == 2):
            cl='red'

        if len(sys.argv)>2 and len(sys.argv[-1]) == 3:
            mode = sys.argv[-1]

        if mode == "xyz":
            ax.scatter(xyzv[:,0], xyzv[:,1], xyzv[:,2], color=cl)
        elif mode == "xyv":
            ax.scatter(xyzv[:,0], xyzv[:,1], xyzv[:,3], color=cl)
        elif mode == "yzv":
            ax.scatter(xyzv[:,1], xyzv[:,2], xyzv[:,3], color=cl)
        elif mode == "xzv":
            ax.scatter(xyzv[:,0], xyzv[:,2], xyzv[:,3], color=cl)
    
    plt.show()
