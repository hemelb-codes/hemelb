import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xyzv = np.loadtxt(sys.argv[1])

    mode = "xyz"

    if len(sys.argv)>2:
        mode = sys.argv[2]

    if mode == "xyz":
        ax.scatter(xyzv[:,0], xyzv[:,1], xyzv[:,2])
    elif mode == "xyv":
        ax.scatter(xyzv[:,0], xyzv[:,1], xyzv[:,3])
    elif mode == "yzv":
        ax.scatter(xyzv[:,1], xyzv[:,2], xyzv[:,3])
    elif mode == "xzv":
        ax.scatter(xyzv[:,0], xyzv[:,2], xyzv[:,3])
    

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()
