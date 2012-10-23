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

    ipFinder = GeometryInletPointFinder(profile.OutputGeometryFile)
    inletPointIndices = ipFinder.GetInletData()

    intersectionFinder = SurfaceIntersectionFinder()
    intersectionFinder.SetFileName(surfaceFile)
    intersectionFinder.SetFileUnitLength(surfaceFileScale)

    for inletId, inlet in enumerate(io for io in profile.Iolets if isinstance(io, Inlet)):
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
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputConnection(cellToPoint.GetOutputPort())
        
        base, gmy = os.path.splitext(profile.OutputGeometryFile)
        writer.SetFileName(base + '.inlet%d.vtp' % inletId)
        writer.Write()

if __name__ == "__main__":
    import sys
    for arg in sys.argv[1:]:
        Run(arg)
    