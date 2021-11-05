# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""
Fabric definitions for HemeLB
Usage:
        With current working directory anywhere inside the HemeLB mercurial checkout, execute
                "fab <machinename> <task>"
        for example:
                "fab hector deploy_cold"

        Do fab -l to get a list of available commands.

        Before use, you MUST copy deploy/machines_user_example.yml as deploy/machines_user.yml and fill in your personal details.
        For smoothest usage, you should also ensure ssh keys are in place so the target machine can access the mercurial repository
        of hemelb.
"""
from templates import *
from machines import *
from fabric.contrib.project import *
from xml.etree import ElementTree
import time
import re
import numpy as np
import yaml
import tempfile
import fab
from os.path import expanduser

# Plotting dependencies are currently in this file.
# TODO: split this off to another child file which only handles plotting?
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import glob

def full_context(cmd):
    return ". ~/.profile && "+cmd

@task
def analyze_stress(results):
    """
    Extracts the Wall shear stress values from the geometry and creates images + ascii files for each extraction.
    """
    #Path to home
    #home = expanduser("~")
    #Path to HemeLB Tools/analysis
    #analysis= home+"/hemelb-dev/hemelb/Tools/analysis/"
    print env.analysis_path
    print template("python $analysis_path/stress.py $results_path/"+results+"/")
    local(template("python $analysis_path/stress.py $results_path/"+results+"/"))

#    os.chdir(env.analysis_path)
#    os.system("python stress.py "env.results_path+results+"/")

@task
def plot_1col_wss_file(filename_filter):
    filenames = glob.glob(filename_filter)
    print filenames
    max_x = np.array([])
    mean_x = np.array([])

    for filename in filenames:
        x = np.loadtxt(filename)
        plt.plot(x)
        if(len(x)>0):
            max_x = np.append(max_x, np.max(x))
            mean_x = np.append(mean_x, np.mean(x))
    plt.ylabel('value')
    plt.show()

    plt.plot(max_x)
    plt.plot(mean_x)
    plt.show()

@task
def plot_4col_file(filename,columns="012",color='blue'):
    """
    Plots an ascii text file with 4 columns as a 3D plot.
    columns="012" indicates that the first three columns are plotted.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    xyzv = np.loadtxt(filename)

    # xyzv = np.where((xyzv[0,:] < -0.208) & (xyzv[2,:] < 0.024))
#    p = xyzv[ (0.023>xyzv[:,0]) & (xyzv[:,0]>0.0225) & (xyzv[:,1]>-0.2109) & (xyzv[:,1]<-0.2107) & (xyzv[:,2]>-0.1548) & (xyzv[:,2]<-0.1544)]

# sample point is: 0.0229437 -0.210787  -0.154438     
#    print p[142]
 
    xyzv = xyzv[1::500]
#    p = p[1::50]
    print len(p)

    if columns == "012":
        ax.scatter(xyzv[:,0], xyzv[:,1], xyzv[:,2], color=color, marker='.')
    elif columns == "013":
        ax.scatter(xyzv[:,0], xyzv[:,1], xyzv[:,3], color=color, marker='.')
    elif columns == "123":
        ax.scatter(xyzv[:,1], xyzv[:,2], xyzv[:,3], color=color, marker='.')
    elif columns == "023":
        ax.scatter(xyzv[:,0], xyzv[:,2], xyzv[:,3], color=color, marker='.')

    ax.scatter(p[:,0], p[:,1], p[:,2], color='red')

    plt.show()


@task
def analyze_velocities(results):

#    local(template("python $analysis_path/allPlanes.py $results_path/"+results+"/"))
#    local(template("python $analysis_path/velocityField.py $results_path/"+results+"/"))
    local(template("python $analysis_path/velocityMag.py $results_path/"+results+"/"))


#    os.system("python allPlanes.py "+results)
#    os.system("python vectorField.py "+results)
#    os.system("python vectorMag.py "+results)
#    os.system("python stress.py "+results)



