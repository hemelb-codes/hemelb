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

def InsideSphere(ijk):
    radius = 10.0
    
    xyz = ijk - 15.5
    
    r2 = np.sum(xyz**2, axis=-1)
    r = np.sqrt(r2)
    # z = r cos theta
    theta = np.arccos(xyz[...,2] / r)
    # tan phi = y/x
    phi = np.arctan2(xyz[..., 1], xyz[..., 0])

    phi_bins = np.mgrid[-np.pi:np.pi:9j]
    theta_bins = np.mgrid[0:np.pi:8j]

    all_points = np.zeros(theta_bins.shape + phi_bins.shape + (3,), dtype=float)
    sintheta = np.sin(theta_bins)
    all_points[..., 0] = sintheta[:, np.newaxis] * np.cos(phi_bins[np.newaxis, :])
    all_points[..., 1] = sintheta[:, np.newaxis] * np.sin(phi_bins[np.newaxis, :])
    all_points[..., 2] = np.cos(theta_bins[:, np.newaxis])
    all_points *= radius
    # a is a point on the plane of the facet
    a = all_points[:-1, :-1]
    b = all_points[1:, :-1]
    c = all_points[0:-1, 1:].copy()
    # Fix the north pole
    c[0, :, :] = all_points[1,1:]
    # n is the normal
    # To avoid any problems, do a cross product
    ac = c-a
    ab = b-a
    n = np.cross(ab, ac)
    n /= np.sqrt(np.sum(n**2, axis=-1))[..., np.newaxis]
    aDOTn = np.sum(a*n,axis=-1)
    
    i_theta = theta_bins.searchsorted(theta) - 1
    i_phi = phi_bins.searchsorted(phi) - 1
    
    xyzDOTn = np.sum(xyz*n[(i_theta, i_phi)], axis=-1)

    return xyzDOTn < aDOTn[(i_theta, i_phi)]
