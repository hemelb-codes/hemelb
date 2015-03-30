from TestResources.sphere import GetSphereNumpy
import numpy as np
from . import TriangleSorter
from .HemeOct import Tree

def mk_trivial():
    # Put the tri points on a square, in the (0,0,0) octant
    points = np.array([(1,1,1),
                       (1,1,2),
                       (1,2,1),
                       (1,2,2)], dtype=float)
    # Define the tris
    triangles = np.array([(0,1,2),
                          (2,1,3)], dtype=int)
    return points, triangles

def test_trivial():
    # 16 cube
    levels = 4
    # put triangles onto the 8 cube level
    tri_level = 3
    points, triangles = mk_trivial()
    
    # Sort!
    tree = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    
    assert tree.levels == levels
    i = tri_level
    for node in tree.IterDepthFirst():
        assert np.all(node.offset == 0)
        assert node.levels == i
        i += 1
        continue
    assert i == levels +1
    
    node = tree.GetNode(tri_level, np.array((0,0,0)))
    ids = node.triIds
    assert np.all(ids == (0,1))


def overlap1d(amin, amax, bmin, bmax):
    return np.logical_and(amin < bmax, amax > bmin)
def overlap3d(amin, amax, bmin, bmax):
    return np.logical_and(
        np.logical_and(
            overlap1d(amin[..., 0], amax[..., 0], bmin[..., 0], bmax[..., 0]),
            overlap1d(amin[..., 1], amax[..., 1], bmin[..., 1], bmax[..., 1])),
        overlap1d(amin[..., 2], amax[..., 2], bmin[..., 2], bmax[..., 2]))
    
def test_sphere():
    levels = 5
    tri_level = 3
    points, triangles, normals = GetSphereNumpy()
    tree = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    
    seen_tris = set()
    for node in tree.IterDepthFirst(tri_level, tri_level):
        n = len(node.triIds)
        node_points = np.zeros((n, 3, 3))
        node_point_ids = triangles[node.triIds]
        for iPt in xrange(3):
            node_points[:,iPt,:] = points[node_point_ids[:,iPt]]
        
        tri_min = node_points.min(axis=1)
        tri_max = node_points.max(axis=1)
        
        nd_min = (node.offset - 1)[np.newaxis,:]
        nd_max = (node.offset + 2**tri_level + 1)[np.newaxis,:]
        
        assert np.all(overlap3d(tri_min, tri_max, nd_min, nd_max))
            
        seen_tris.update(node.triIds)
    allids = sorted(seen_tris)
    assert np.all(np.equal(allids, np.arange(len(triangles))))
    