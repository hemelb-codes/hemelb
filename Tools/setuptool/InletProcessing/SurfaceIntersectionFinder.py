# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from vtk import vtkSTLReader, vtkCutter, vtkPlane, vtkPolyDataConnectivityFilter, vtkTransformPolyDataFilter, vtkTransform

class SurfaceIntersectionFinder(object):
    """Compute the edge of the surface at the inlet.
    """
    def __init__(self, surfaceFile, fileUnitLength):
        """surfaceFile - the path of the STL file containing the surface
        fileUnitLength - the length, in metres, of 1 in the STL file
        """
        self.reader = vtkSTLReader()
        self.reader.SetFileName(surfaceFile)
        
        self.cutter = vtkCutter()
        self.cutter.SetInputConnection(self.reader.GetOutputPort())
        
        self.regionPicker = vtkPolyDataConnectivityFilter()
        self.regionPicker.SetInputConnection(self.cutter.GetOutputPort())
        self.regionPicker.SetExtractionModeToClosestPointRegion()
        
        self.scaler = vtkTransformPolyDataFilter()
        self.scaler.SetInputConnection(self.regionPicker.GetOutputPort())
        trans = vtkTransform()
        self.fileUnitLength = fileUnitLength
        trans.Scale(fileUnitLength,fileUnitLength,fileUnitLength)
        self.scaler.SetTransform(trans)
        
    def IntersectWithInlet(self, inlet):
        """Return the intersection of the surface with the supplied inlet. 
        """
        plane = vtkPlane()
        r = inlet.Centre
        plane.SetOrigin(r.x, r.y, r.z)
        n = inlet.Normal
        plane.SetNormal(n.x, n.y, n.z)
        self.cutter.SetCutFunction(plane)
        
        self.regionPicker.SetClosestPoint(r.x, r.y, r.z)
        
        self.scaler.Update()
        r.x *= self.fileUnitLength
        r.y *= self.fileUnitLength
        r.z *= self.fileUnitLength
        
        return self.scaler.GetOutput()
    