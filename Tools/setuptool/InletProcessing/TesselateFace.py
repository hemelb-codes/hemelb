import numpy as np
from vtk import vtkDelaunay2D, vtkTransformPolyDataFilter, vtkPoints, \
    vtkPolyData, vtkTransform

class Tesselator(object):
    """Object to transform the iolet to the X-Y plane and tesselate the 
    interior with a triangular mesh.
    """

    def __init__(self, inlet, intersectionPD, positionsNP):
        """Create the object that will tesselate the face.
        
        inlet -- the HemeLbSetupTool.Model.Iolets.Iolet object
        
        intersectionPD -- vtkPolyData representing the intersection between the 
            iolet and the surface
            
        positionsNP -- numpy array of the positions of the LB sites that are
            adjacent to the iolet
        """
        self.trans = self._CreateTransformToInletCoords(inlet)
        mergedPD = self._MergePoints(intersectionPD, positionsNP)
        
        tf = vtkTransformPolyDataFilter()
        tf.SetTransform(self.trans)
        tf.SetInput(mergedPD)
        
        self.tesselator = vtkDelaunay2D()
        self.tesselator.SetSource(intersectionPD)
        self.tesselator.SetInputConnection(tf.GetOutputPort())
        return
    
    def __call__(self):
        self.tesselator.Update()
        return self.tesselator.GetOutput()
    
    @staticmethod
    def _MergePoints(intersectionPD, positionsNP):
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
    
    @staticmethod
    def _CreateTransformToInletCoords(inlet):
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
