# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import os.path
import numpy as np
import vtk

from HemeLbSetupTool.Model.Profile import Profile
from HemeLbSetupTool.Model.Iolets import Inlet

from InletPointFinder import GeometryInletPointFinder
from SurfaceIntersectionFinder import SurfaceIntersectionFinder
from TesselateFace import Tesselator
from PoiseuilleSolver import PoiseuilleSolver

def Run(profileFile):
    """Process all the Inlets specified by profileFile, solving the Navier-
    Stokes equations for steady flow down an infinite channel with the cross
    section of the inlet. Writes the point ID and the fluid velocity at that 
    point, the velocity being such that the integral of the flow across the
    channel is unity.
    
    Currently output is written as vtkPolyData to files like
    "$GEOMETRYFILEBASENAME.inlet$ID.vtp"
    """
    # Load the profile
    profile = Profile()
    profile.LoadFromFile(profileFile)
    
    surfaceFile = profile.StlFile
    surfaceFileScale = profile.StlFileUnit.SizeInMetres

    print profile.OutputXmlFile
    ipFinder = GeometryInletPointFinder(profile.OutputXmlFile)
    inletPointIndices = ipFinder.GetInletData()
    if len(inletPointIndices) == 0:
        print "WARNING: no inletPointIndices found. This may mean that HemeLb is run with no inlets inside the simulation domain (e.g., the inlets defined in the xml file reside outside of the physical geometry)."

    intersectionFinder = SurfaceIntersectionFinder()
    intersectionFinder.SetFileName(surfaceFile)
    intersectionFinder.SetFileUnitLength(surfaceFileScale)

    for inletId, inlet in enumerate(io for io in profile.Iolets if isinstance(io, Inlet)):
        print inletId, len(profile.Iolets), len(inletPointIndices)
        inletPointPD = inletPointIndices[inletId]
        intersectionFinder.SetIolet(inlet)
        tesselator = Tesselator()
        tesselator.SetInlet(inlet)
        tesselator.SetEdgeConnection(intersectionFinder.GetOutputPort())
        tesselator.SetSitesConnection(inletPointPD.GetProducerPort())
        
        solver = PoiseuilleSolver()
        solver.SetInputConnection(tesselator.GetOutputPort())
        
        cellToPoint = vtk.vtkCellDataToPointData()
        cellToPoint.SetInputConnection(solver.GetOutputPort())
        cellToPoint.Update()
        #writer = vtk.vtkXMLPolyDataWriter()
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileTypeToASCII()
        writer.SetInputConnection(cellToPoint.GetOutputPort())
       
        # print type(cellToPoint) #cellToPoint is of type vtkobject
        # print dir(cellToPoint)
 
        base, gmy = os.path.splitext(profile.OutputGeometryFile)
        writer.SetFileName(base + '.inlet%d.txt' % inletId)
        writer.Write()
        return base + '.inlet%d.txt' % inletId

def CreateInletWeightsFile(fname):
    """ This function creates a weighting file for the inlet. It 
    takes a textual output file and reorganizes this data into a 
    simpler and more compact form. The format is as follows:

    <lattice x coord> <lattice y coord> <lattice z coord> <weight>
    """

    RecordVelocities = False

    vels = np.zeros((1))
    offset_counter = 0

    f = open(fname)
    for line in f:

        if RecordVelocities == True:
            vels_line = line.strip().split()
            for v in vels_line:
              vels[offset_counter] = float(v)
              offset_counter += 1

        # Line preceding all the velocity values
        if "LOOKUP_TABLE default" in line:
            RecordVelocities = True
        
        if "POINT_DATA" in line:
            p = line.strip().split()
            print "Number of elements is: ", p[1]
            vels = np.zeros((int(p[1])))

    print vels

if __name__ == "__main__":
    import sys
    import time
    for arg in sys.argv[1:]:
        start = time.time()
        txtOutputFileName = Run(arg)
        CreateInletWeightsFile(txtOutputFileName)
        end = time.time()
        print end-start," seconds elapsed."
    
