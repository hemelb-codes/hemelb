from vtk import vtkSTLReader, vtkCutter, vtkPlane, \
    vtkPolyDataConnectivityFilter, vtkTransformPolyDataFilter, vtkTransform

class SurfaceIntersectionFinder(object):
    """Contains a VTK pipeline to find the intersection of a Iolet plane with
    the contents of an STL file, and return this data in SI units.
    
    Pipeline is:
    
    FileUnitLength ----------------> vtkTransform ------------------------
                                                                          \
    FileName -> vtkSTLReader -> vtkCutter -> vtkPolyDataConnectivityFilter -> vtkTransformPolyDataFilter
                            /            /
    Iolet -----> vtkPlane -            /
            \-------- Centre ---------
    
    """
    
    def __init__(self):
        """Setup the pipeline.
        """
        self.reader = vtkSTLReader()
        
        self.cutter = vtkCutter()
        self.cutter.SetInputConnection(self.reader.GetOutputPort())
        
        self.regionPicker = vtkPolyDataConnectivityFilter()
        self.regionPicker.SetInputConnection(self.cutter.GetOutputPort())
        self.regionPicker.SetExtractionModeToClosestPointRegion()
        
        self.scaler = vtkTransformPolyDataFilter()
        self.scaler.SetInputConnection(self.regionPicker.GetOutputPort())
        
        self.GetFileName = self.reader.GetFileName
        self.SetFileName = self.reader.SetFileName
        
        self.GetOutputPort = self.scaler.GetOutputPort
        self.GetOutput = self.scaler.GetOutput
        self.Update = self.scaler.Update
        
        self.iolet = None
        self.fileUnitLength = None
        
        return
        
    def SetFileUnitLength(self, fileUnitLength):
        """Set the length, in metres, of 1 in the STL file.
        """
        self.fileUnitLength = fileUnitLength
        trans = vtkTransform()
        trans.Scale(fileUnitLength,fileUnitLength,fileUnitLength)
        self.scaler.SetTransform(trans)
        return
    def GetFileUnitLength(self):
        return self.fileUnitLength
    
    
    def SetIolet(self, iolet):
        self.iolet = iolet
        
        plane = vtkPlane()
        r = iolet.Centre
        plane.SetOrigin(r.x, r.y, r.z)
        n = iolet.Normal
        plane.SetNormal(n.x, n.y, n.z)
        self.cutter.SetCutFunction(plane)
        
        self.regionPicker.SetClosestPoint(r.x, r.y, r.z)
        return
    
    def GetIolet(self):
        return self.iolet
    pass