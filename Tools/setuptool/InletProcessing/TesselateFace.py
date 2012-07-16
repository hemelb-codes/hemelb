# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import numpy as np
from vtk import vtkDelaunay2D, vtkTransformPolyDataFilter

def TransformAndTesselate(inlet, intersectionPD, positionsNP):
    trans = CreateTransformToInletCoords(inlet)
    mergedPD = MergePoints(intersectionPD, positionsNP)
    
    tf = vtkTransformPolyDataFilter()
    tf.SetTransform(trans)
    tf.SetInput(mergedPD)
    
    tesselator = vtkDelaunay2D()
    tesselator.SetSource(intersectionPD)
    tesselator.SetInputConnection(tf.GetOutputPort())
    
    tesselator.Update()
    return tesselator.GetOutput()

from vtk import vtkPoints, vtkPolyData
def MergePoints(intersectionPD, positionsNP):
    """Create a vtkPolyData containing all the points, with the ones defining 
    the perimeter first and in the same order as the input.
    """
    nI = intersectionPD.GetNumberOfPoints()
    nP = len(positionsNP)
    
    outPD = vtkPolyData()
    outPoints = vtkPoints()
    outPoints.DeepCopy(intersectionPD.GetPoints())
    
    outPoints.SetNumberOfPoints(nI + nP)
    for iPos, iOut in enumerate(xrange(nI, nI + nP)):
        outPoints.SetPoint(iOut, *positionsNP[iPos])
        continue
    
    outPD.SetPoints(outPoints)
    return outPD

from vtk import vtkTransform
def CreateTransformToInletCoords(inlet):
    """Create a vtkTransform which will rotate the coordinates such that
    inlet.Normal is oriented in the z-direction and origin is inlet.Centre
    """
    n = inlet.Normal
    n = np.array([n.x, n.y, n.z])
    
    z = np.array([0., 0., 1.])
    
    theta = np.arccos(np.dot(n,z))
    axis  = np.cross(n, z)
    
    norm = np.sqrt(np.dot(axis, axis))
    axis /= norm
    
    trans = vtkTransform()
    trans.Scale(1., 1., 0.)
    trans.RotateWXYZ(180. * theta / np.pi, axis[0], axis[1], axis[2])
    
    r = inlet.Centre
    trans.Translate(r.x, r.y, r.z)
    
    return trans
