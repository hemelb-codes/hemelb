# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from HemeLbSetupTool.Model.Profile import Profile
from HemeLbSetupTool.Model.Iolets import Inlet
from InletPointFinder import InletPointFinder
from SurfaceIntersectionFinder import SurfaceIntersectionFinder

def IterInletsData(profileFile):
    """Iterator over the INLETS, returning the inlet lattice point 3D indices 
    and the line defining the edge of the domain.
    """
    profile = Profile()
    profile.LoadFromFile(profileFile)
    
    geometryFile = profile.OutputGeometryFile
    surfaceFile = profile.StlFile
    surfaceFileScale = profile.StlFileUnit.SizeInMetres
    
    inlets = filter(lambda io: isinstance(io, Inlet), profile.Iolets)
    ipFinder = InletPointFinder(geometryFile)
    inletPointIndices = ipFinder.Run()
    
    intersectionFinder = SurfaceIntersectionFinder(surfaceFile, surfaceFileScale)
    
    for i, inlet in enumerate(inlets):
        yield inlet, inletPointIndices[i], intersectionFinder.IntersectWithInlet(inlet)
    
