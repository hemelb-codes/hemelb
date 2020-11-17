# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

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
        
        #solver = PoiseuilleSolver()
        #solver.SetInputConnection(tesselator.GetOutputPort())
        
        cellToPoint = vtk.vtkCellDataToPointData()
        cellToPoint.SetInputConnection(tesselator.GetOutputPort())
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

def CreateEmptyInletWeightsFile(fname_in, fname_out):
    """ This function creates an empty weighting file for the inlet. It 
    takes a textual output file and reorganizes this data into a 
    simpler and more compact form. The format is as follows:

    <lattice x coord> <lattice y coord> <lattice z coord> 0.0
    """

    RecordCoordinates = False

    coords = np.zeros ((3))
    coords_offset_counter = 0
 
    # READING STEP
    f = open(fname_in)
    for line in f:

        if RecordCoordinates == True:
            """
            This branch extracts all the coordinate values from the VTK version 3 format file.
            """
            coords_line = line.strip().split()
            for c in coords_line:
                coords[coords_offset_counter] = int(c)
                coords_offset_counter += 1
                if coords_offset_counter == len(coords):
                    RecordCoordinates = False
                    break

        elif "POINT_DATA" in line:
            """
            Line indicating the number of velocity values.
            Parsing condition looks silly because the VTK version 3 format is generic, but 
            not very descriptive.
            """
            p = line.strip().split()
            print "Number of elements is: ", int(p[1])

            coords = np.zeros((int(p[1])*3))

        elif "PointIds" in line:
            """
            TODO: detect the coordinate values.
            """
            RecordCoordinates = True

    coords = coords.reshape(len(coords)/3, 3)

    # WRITING STEP
    f = open("%s" % (fname_out),'w')

    for i in xrange(0, len(coords)):
        # take out faulty data points at the origin or with -1,-1,-1 coordinates.
        if sum(coords[i])>0.5:
            f.write("%i %i %i 0.0\n" % (int(coords[i][0]), int(coords[i][1]), int(coords[i][2])))


def CreateInletWeightsFile(fname):
    """ This function creates a weighting file for the inlet. It 
    takes a textual output file and reorganizes this data into a 
    simpler and more compact form. The format is as follows:

    <lattice x coord> <lattice y coord> <lattice z coord> <weight>
    """

    RecordVelocities = False
    RecordCoordinates = False

    vels = np.zeros((1)) # NP Array which stores all the original velocities (they add up to 1 supposedly)
    coords = np.zeros ((3))
    vels_offset_counter = 0
    coords_offset_counter = 0
    
 
    f = open(fname)
    for line in f:

        if RecordVelocities == True:
            """
            This branch extracts all the velocity values from the VTK version 3 format file.
            """
            vels_line = line.strip().split()
            for v in vels_line:
                print "offset = ", vels_offset_counter, v
                vels[vels_offset_counter] = float(v)
                vels_offset_counter += 1
                if vels_offset_counter == len(vels):
                    RecordVelocities = False
                    break

        elif RecordCoordinates == True:
            """
            This branch extracts all the coordinate values from the VTK version 3 format file.
            """
            coords_line = line.strip().split()
            for c in coords_line:
                coords[coords_offset_counter] = int(c)
                coords_offset_counter += 1
                if coords_offset_counter == len(coords):
                    RecordCoordinates = False
                    break

        elif "LOOKUP_TABLE default" in line:
            """
            Line preceding all the velocity values
            """
            RecordVelocities = True
        
        elif "POINT_DATA" in line:
            """
            Line indicating the number of velocity values.
            Parsing condition looks silly because the VTK version 3 format is generic, but 
            not very descriptive.
            """
            p = line.strip().split()
            print "Number of elements is: ", int(p[1])

            vels   = np.zeros((int(p[1])))
            coords = np.zeros((int(p[1])*3))

        elif "PointIds" in line:
            """
            TODO: detect the coordinate values.
            """
            RecordCoordinates = True

    """
    We renormalize the velocities, so that the maximum velocity is equal to 1.
    From now on: we use these values now as _velocity weights_.
    """
    print len(vels)
    coords = coords.reshape(len(coords)/3, 3)

    vels = vels * (1.0 / vels.max()) 

    f = open("%s.weights.txt" % (fname[:-4]),'w')
    for i in xrange(0, len(vels)):
        # take out faulty data points at the origin or with -1,-1,-1 coordinates.
        if sum(coords[i])>0.5:
            f.write("%i %i %i %f\n" % (int(coords[i][0]), int(coords[i][1]), int(coords[i][2]), vels[i]))

if __name__ == "__main__":
    import sys
    import time
    start = time.time()

    mode = "emptyweights"

    if mode == "full":
        """ Take a profile file, gmy + input, and create a weights file """
        txtOutputFileName = Run(sys.argv[1])
        CreateInletWeightsFile(txtOutputFileName)
    elif mode == "weight_only":
        """ Take an intermediate VTK version 3 ASCII format file and convert it to a smaller inlets weights file. """
        CreateInletWeightsFile(sys.argv[1])
    if mode == "emptyweights":
        """ Take a profile file, gmy + input, and create a weights file """
        if len(sys.argv) < 3:
            print "Usage: python <script_name> <input pr2 file name> <output empty weights file name>"
            sys.exit()
        txtOutputFileName = Run(sys.argv[1])
        CreateEmptyInletWeightsFile(txtOutputFileName, sys.argv[2])
    end = time.time()
    print end-start," seconds elapsed."
    

