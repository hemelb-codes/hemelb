import vtk.util.numpy_support
from hemeTools.utils.cd import cd
import os.path
import numpy as np

_sphere = None
def GetSphere():
    """The sphere has radius 10 and is centred on 15.5. It has normals.
    It requires a box with 5 levels.
    """ 
    global  _sphere
    if _sphere is None:
        with cd(os.path.dirname(__file__)):
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName('sphere.vtp')
            reader.Update()
            _sphere = reader.GetOutput()
    return _sphere

_np_sphere = None
def GetSphereNumpy():
    """The sphere data in Numpy arrays.
    """
    global _np_sphere
    if _np_sphere is None:
        surf = GetSphere()
        tridata = vtk.util.numpy_support.vtk_to_numpy(surf.GetPolys().GetData())
        tridata.shape = (surf.GetNumberOfPolys(), 4)
        assert np.all(tridata[:,0]==3), "Non triangle in surface!"
        triangles = tridata[:, 1:]
        points = vtk.util.numpy_support.vtk_to_numpy(surf.GetPoints().GetData())
        normals = surf.GetCellData().GetNormals()
        normals = vtk.util.numpy_support.vtk_to_numpy(normals)
        _np_sphere = (points, triangles, normals)
    return _np_sphere

def Radius(ijk):
    xyz = ijk - [15.5, 15.5, 15.5]
    r2 = np.sum(xyz**2, axis=-1)
    return np.sqrt(r2)

def InsidePerfectSphere(ijk):
    """Return true if the point is inside the (perfect) sphere.
    """
    return Radius(ijk) < 10.0