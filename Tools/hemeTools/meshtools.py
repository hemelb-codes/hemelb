import os.path
import vtk
import numpy as N
from . import cache

# These coords are in the HemeLB basis I.e., same axes, but origin at
# the centre of mass of the triangles (assuming equal tri mass).
# Hence need to:
# - compute centroid of each tri, then average centroids
# - add this to each *let centre.

def readSTLSurfaceFile(stlFile):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(stlFile)
    reader.Update()
    return reader.GetOutput()

@cache.processesFile
def processMesh(stlFile):
    surface = readSTLSurfaceFile(stlFile)
    
    tri_ctr_x = N.zeros(3, dtype=N.float)
    
    min_x = N.zeros(3)
    max_x = N.zeros(3)
    ctr_x = N.zeros(3)

    nTri = surface.GetNumberOfCells()
    for i in xrange(nTri):
        tri = surface.GetCell(i).GetPoints()
        tri_ctr_x[:] = 0.
        # for each vertex of tri i
        for m in xrange(3):
            tri_ctr_x += tri.GetPoint(m)
            continue
        continue
    
    for m in range(3):
        p = tri_ctr_x + 1.01 * (tri.GetPoint(m) - tri_ctr_x)
        tri.SetPoint(m, *p)
        
        min_x = N.fmin(min_x, p)
        min_x = N.fmin(min_x, p)
        
        ctr_x += p
        continue
    
    ctr_x /= 3.*nTri
    
    voxel_size = 0.
    for n in xrange(nTri):
        
        tri = surface.GetCell(n).GetPoints()
        for m in range(3):
            tri.SetPoint(m,
                         *(tri.GetPoint(m) - ctr_x))
            continue
        
        dx = N.array(tri.GetPoint(1)) - N.array(tri.GetPoint(0))
        voxel_size += N.sum(dx**2)
        
        dx = N.array(tri.GetPoint(2)) - N.array(tri.GetPoint(0))
        voxel_size += N.sum(dx**2)
        
        dx = N.array(tri.GetPoint(2)) - N.array(tri.GetPoint(1))
        voxel_size += N.sum(dx**2)
        continue
    
    voxel_size /= nTri * 3
    
    voxels = N.ceil(1.5 * (max_x - min_x)/voxel_size).astype(N.int)
    dim = voxel_size * voxels
    half_dim = 0.5 * dim
    
    out = {}
    outvars = ['ctr_x', 'voxel_size', 'voxels', 'dim', 'min_x', 'max_x']
    for var in outvars:
        out[var] = locals()[var]
        continue
    return out

def meshCentre(stlFile):
    out = processMesh(stlFile)
    return out['ctr_x']
