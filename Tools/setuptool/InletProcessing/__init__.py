import numpy as np
from GatherInletData import IterInletsData
from TesselateFace import Tesselator
from PoiseuilleSolver import PoiseuilleSolver
import vtk
from vtk.util import numpy_support as convert

import pdb
def DbgWrite(polyDataOrAlgorithm, filename):
    w = vtk.vtkXMLPolyDataWriter()
    if isinstance(polyDataOrAlgorithm, vtk.vtkAlgorithm):
        w.SetInputConnection(polyDataOrAlgorithm.GetOutputPort())
    elif isinstance(polyDataOrAlgorithm, vtk.vtkPolyData):
        w.SetInput(polyDataOrAlgorithm)
    else:
        raise ValueError('cannae cope with a "%s"' % str(type(polyDataOrAlgorithm)))
    
    w.SetFileName(filename)
    w.Write()

def Run(profileFile):
    for inletId, (inlet, (ids, positions), intersection) in enumerate(IterInletsData(profileFile)):
        tesselator = Tesselator(inlet, intersection, positions)
        tesselated = tesselator()
        DbgWrite(tesselated, 'tess.vtp')
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
        
