# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from math import sqrt
from vtk import vtkProgrammableFilter, vtkPoints, vtkFloatArray, vtkPolyData

class AverageSideLengthCalculator(vtkProgrammableFilter):
    """Given an input vtkPolyData object, output a vtkPolyData with
    one point whose associated scalar is the average side length of
    the triangle contained within the input.

    Included is a convenience method GetOutputValue which will return
    a Python float of the answer.
    
    """
    def __init__(self):
        self.SetExecuteMethod(self.Execute)
        
    def Execute(self, *args):
        polydata = self.GetPolyDataInput()
        nTris = polydata.GetNumberOfCells()
        
        totalPerim = 0.
        for i in xrange(nTris):
            tri = polydata.GetCell(i)
            perim = 0.
            for j in range(3):
                perim += sqrt(tri.GetEdge(j).GetLength2())
                continue
            totalPerim += perim
            continue
        
        aveSide = totalPerim / (3*nTris)

        p = vtkPoints()
        p.InsertPoint(0, 0.,0.,0.)
        
        val = vtkFloatArray()
        val.InsertValue(0, aveSide)
        
        out = vtkPolyData()
        out.SetPoints(p)
        out.GetPointData().SetScalars(val)
        
        self.GetPolyDataOutput().ShallowCopy(out)
        return

    def GetOutputValue(self):
        self.Update()
        vals = self.GetOutput().GetPointData().GetScalars()
        if vals is None:
            return None
        return vals.GetValue(0)
    
    pass
