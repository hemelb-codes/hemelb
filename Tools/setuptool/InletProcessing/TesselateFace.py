# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import numpy as np
from vtk import vtkDelaunay2D, vtkTransformPolyDataFilter, vtkTransform, \
    vtkProgrammableFilter, vtkAppendPolyData, vtkIntArray


class Tesselator(object):
    """Object to transform the iolet sites and edge to the X-Y plane and 
    tesselate the interior with a triangular mesh.
    
    Pipeline is:
    
    Edge ---------------------------------------------------
          \                                                 \
    Sites--> vtkAppendPolyData -> vtkTransformPolyDataFilter -> vtkDelaunay2D
                              /    
                            /
    Iolet -> vtkTransform-/
    """

    def __init__(self):
        self.dummyAdder = DummyPointIdsAdder()
        
        self.merger = vtkAppendPolyData()
        self.merger.AddInputConnection(self.dummyAdder.GetOutputPort())
        
        self.transformer = vtkTransformPolyDataFilter()
        self.transformer.SetInputConnection(self.merger.GetOutputPort())
        
        self.tesselator = vtkDelaunay2D()
        self.tesselator.SetInputConnection(self.transformer.GetOutputPort())
        
        self.GetOutputPort = self.tesselator.GetOutputPort
        self.GetOutput = self.tesselator.GetOutput
        self.Update = self.tesselator.Update
        
        return
    
    def SetInlet(self, inlet):
        """Create a vtkTransform which will rotate the coordinates such that
        inlet.Normal is oriented in the z-direction and origin is inlet.Centre
        """
        n = inlet.Normal
        n = np.array([n.x, n.y, n.z])
        
        minInd = n.argsort()[0]
        axis = np.zeros(3)
        axis[minInd] = 1
        
        transmat = np.eye(4)
        transmat[0, 0:3] = np.cross(n, axis)
        transmat[0, 0:3] /= np.sqrt(np.sum(transmat[0, 0:3]**2))
        transmat[1, 0:3] = np.cross(n, transmat[0, 0:3])
        transmat[2, 0:3] = n
        
        trans = vtkTransform()
        trans.Scale(1., 1., 0.)
        trans.Concatenate(transmat.flatten())
        r = inlet.Centre
        trans.Translate(r.x, r.y, r.z)
        
        self.trans = trans
        self.transformer.SetTransform(self.trans)
        return
    
    def SetEdgeConnection(self, edgePort):
        self.dummyAdder.SetInputConnection(edgePort)
        self.tesselator.SetSourceConnection(edgePort)
        return
    
    def SetSitesConnection(self, sitesPort):
        self.merger.AddInputConnection(sitesPort)
        return
    pass

class DummyPointIdsAdder(vtkProgrammableFilter):
    """Filter to add PointData with an array of PointIds to the edge 
    vtkPolyData. This is needed because vtkAppendPolyData only passes through
    fields that are present in ALL input PolyData.
    
    Use (-1,-1,-1) as the dummy value (since (0,0,0) is the lowest expected
    in a .gmy file).
    """
    def __init__(self):
        self.SetExecuteMethod(self._Execute)
        return
    
    def _Execute(self):
        input = self.GetPolyDataInput()
        output = self.GetPolyDataOutput()
        output.ShallowCopy(input)
        
        empty = vtkIntArray()
        empty.SetNumberOfComponents(3)
        empty.SetNumberOfTuples(input.GetNumberOfPoints())
        for i in xrange(3):
            empty.FillComponent(i, -1)
        empty.SetName("PointIds")
            
        opd = output.GetPointData()
        opd.AddArray(empty)
        return
    
    pass
