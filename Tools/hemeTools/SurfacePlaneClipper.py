import vtk
import numpy as N
from .Plane import Plane

def planeToVTK(pp):
    """Create an equivalent vtkPlane object."""
    vp = vtk.vtkPlane()
    if pp.normal is not None:
        vp.SetNormal(pp.normal)
    if pp.centre is not None:
        vp.SetOrigin(pp.centre)
    return vp


def normalDot(p1, p2):
    """Dot product between two planes, either vtkPlane or
    PlaneClipper.Plane instances.
    
    """
    n1 = p1.GetNormal()
    n2 = p2.GetNormal()
    return (n1[0]*n2[0] +
            n1[1]*n2[1] +
            n1[2]*n2[2])

class SurfacePlaneClipper(object):
    """A VMTK-like module, for clipping a surface against one or more
    planes. The result will be stored in the instance variable
    Output."""
    
    def __init__(self, Surface=None, Planes=None):
        """surface -- The vtkPolyData object we want to clip.
        
        planes -- A container of the planes against which to clip,
        must be vtkPlane or PlaneClipper.Plane instances.

        Both can be omitted here, but must be set on the instance
        before Execute() is called.
        """
        self.Surface = Surface
        self.Planes = Planes
        self.Clipper = None
        self.Clipped = None
        self.Output = None
        return 
    
    def clip(self, plane):
        """Clip the surface against the supplied plane, using a
        vtkClipPolyData filter."""
        
        if self.Clipper is None:
            self.Clipper = vtk.vtkClipPolyData()
            self.Clipper.SetInput(self.Clipped)
            pass
        
        # Get the axis-aligned cuboid bounding the surface
        bounds = self.Surface.GetBounds()
        # (xmin,xmax, ymin,ymax, zmin,zmax)
        planes = vtk.vtkPlanes()
        planes.SetBounds(*bounds)
        # Iterate over BB, find plane with closest normal, replace
        # with our plane.
        maxDot = -1; maxI = -1
        for i in range(planes.GetNumberOfPlanes()):
            bound = planes.GetPlane(i)
            dot = normalDot(plane, bound)
            if dot> maxDot:
                maxDot = dot
                maxI = i
                pass
            continue
        
        planes.GetNormals().SetTuple3(maxI, *plane.GetNormal())
        planes.GetPoints().SetPoint(maxI, *plane.GetOrigin())
        
        self.Clipper.SetClipFunction(planes)
        self.Clipper.GenerateClippedOutputOn()
        self.Clipper.InsideOutOn()
        self.Clipper.Update()
        self.Clipped.DeepCopy(self.Clipper.GetClippedOutput())
        self.Clipped.Update()
        return    
    
    
    def GetTargetSurface(self):
        """Get the connected part of the clipped surface that lies
        closest to the first clipping plane's origin, using
        vtkPolyDataConnectivityFilter.
        
        """
        filter = vtk.vtkPolyDataConnectivityFilter()
        filter.SetExtractionModeToClosestPointRegion()
        filter.SetClosestPoint(self.Planes[0].GetOrigin())
        filter.SetInput(self.Clipped)
        filter.Update()
        self.Output = filter.GetOutput()
        return
        
    def Execute(self):
        """Run the proceedure:
        
        For each plane in Planes, clip Surface against it.  Then
        extract the connected surface that lies closest to the first
        specified plane origin.

        Result is stored in Output.
        
        """
        self.Clipped = getattr(vtk, self.Surface.GetClassName())()
        self.Clipped.DeepCopy(self.Surface)
        
        for p in self.Planes:
            if isinstance(p, Plane):
                p = planeToVTK(p)
                pass
            self.clip(p)
            continue
        
        self.GetTargetSurface()
        return
    
    pass



def clip(infile, planes, outfile):
    import vmtk
    
    reader = vmtk.vmtksurfacereader.vmtkSurfaceReader()
    reader.InputFileName = infile
    reader.Execute()
    
#     planes = [Plane(normal=p.normal,
#                     centre=p.centre+ctr,
#                     name=p.name) for p in planes]
    
    clipper = SurfacePlaneClipper()
    clipper.Surface = reader.Output
    clipper.Planes = planes
    clipper.Execute()
    
    writer = vmtk.vmtksurfacewriter.vmtkSurfaceWriter()
    writer.Surface = clipper.Output
    writer.OutputFileName = outfile
    writer.Execute()
    
    for p in planes:
        print p
    
        
