from vtk import *
import numpy as np


class TriangulatedCylinderSource(vtkProgrammableSource):
    def __init__(self):
        self.Center = (0,0,0)
        self.Height = 1.
        
        self.Centerline = vtkPolyData()
        self.Tuber = vtkTubeFilter()
        self.Triangulator = vtkTriangleFilter()
        self.Triangulator.SetInputConnection(self.Tuber.GetOutputPort())
        self.SetExecuteMethod(self.Execute)
        return
    
    def GetResolution(self):
        return self.Tuber.GetNumberOfSides()
    def SetResolution(self, n):
        return self.Tuber.SetNumberOfSides(n)
    
    def GetCenter(self):
        return self.Center
    def SetCenter(self, c):
        self.Center = c
        return
    
    def GetRadius(self):
        return self.Tuber.GetRadius()
    def SetRadius(self, r):
        return self.Tuber.SetRadius(r)
        
    def GetCapping(self):
        return self.Tuber.GetCapping()
    def SetCapping(self, c):
        return self.Tuber.SetCapping(c)
    def CappingOn(self):
        return self.Tuber.CappingOn()
    def CappingOff(self):
        return self.Tuber.CappingOff()

    def GetHeight(self):
        return self.Height
    def SetHeight(self, h):
        self.Height = h
        return

    def _GetDx(self):
        return 2 * self.GetRadius() * np.sin(np.pi / self.GetResolution())
    def _GetNz(self):
        targetDz = self._GetDx()
        return int(np.round(self.GetHeight()/ targetDz))
    def _GetDz(self):
        return self.GetHeight() / self._GetNz()
    
    def Execute(self, *args):
        self.Centerline.PrepareForNewData()

        nSegments = self._GetNz()
        h = self.GetHeight()
        dz = self._GetDz()
        c = self.GetCenter()
        
        points = vtkPoints()
        points.SetNumberOfPoints(nSegments + 1)

        line = vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(nSegments + 1)

        for i in xrange(nSegments + 1):
            points.SetPoint(i, c[0], c[1], c[2] - 0.5*h + i*dz)
            line.GetPointIds().SetId(i, i)
            continue

        self.Centerline.Allocate(1,1)
        self.Centerline.SetPoints(points)
        self.Centerline.InsertNextCell(line.GetCellType(),
                                       line.GetPointIds())
        
        self.Tuber.SetInput(self.Centerline)
        self.Triangulator.Update()

        out = self.GetPolyDataOutput()
        out.ShallowCopy(self.Triangulator.GetOutput())
        return
    
    pass

