# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import numpy as np
from GatherInletData import IterInletsData
from TesselateFace import TransformAndTesselate
from PoiseuilleSolver import PoiseuilleSolver
import vtk
from vtk.util import numpy_support as convert

import pdb

def Run(profileFile):
    for inletId, (inlet, (ids, positions), intersection) in enumerate(IterInletsData(profileFile)):
        tesselated = TransformAndTesselate(inlet, intersection, positions)
        p = tesselated.GetPoints()
        nPoints = tesselated.GetNumberOfPoints()
        nCells = tesselated.GetNumberOfPolys()
        
#        r = [0., 0., 0.]
#        for i in xrange(nPoints):
#            p.GetPoint(i, r)
#            r[2] = 0.
##            rho.SetTuple1(i, np.sqrt(r[0]**2 + r[1]**2))
#            p.SetPoint(i, r)
#            continue
        
#        tesselated.GetPointData().SetScalars(rho)
        
        pdb.set_trace()
        solver = PoiseuilleSolver(tesselated)
        uField = solver.Solve()
        
        u = convert.numpy_to_vtk(uField.getValue())
        tesselated.GetCellData().SetScalars(u)
        
        w = vtk.vtkXMLPolyDataWriter()
        w.SetInput(tesselated)
        w.SetFileName('tess%d.vtp' % inletId)
        w.Write()
        
