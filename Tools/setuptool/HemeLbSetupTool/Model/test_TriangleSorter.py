import numpy as np
from TestResources.sphere import GetSphereNumpy
from . import TriangleSorter
from .HemeOct import Tree
from nose.tools import set_trace
def trees_with_triIds_equal(t1, t2, tri_level):
    if t1 != t2:
        return False
    # If here, they're structurally the same, check triIds
    for n1 in t1.IterDepthFirst(tri_level, tri_level):
        n2 = t2.GetNode(n1.levels, n1.offset)
        if not np.all(n1.triIds == n2.triIds):
            return False
        continue
    # None of the nodes differ so we are true!
    return True
        
def mk_trivial():
    # Put the tri points on a square, in the (0,0,0) octant
    points = np.array([(1.2, 1.2, 1.2),
                       (1.2, 1.2, 2.2),
                       (1.2, 2.2, 1.2),
                       (1.2, 2.2, 2.2)], dtype=float)
    # Define the tris
    triangles = np.array([(0,1,2),
                          (2,1,3)], dtype=int)
    normals = np.array([[-1, 0, 0],
                        [-1, 0, 0]], dtype=float)
    return points, triangles, normals

def mk_trivial2():

    points = np.array([(1.2, 1.2, 1.2),
                       (1.2, 1.2, 2.2),
                       (1.2, 2.2, 1.2),
                       (1.2, 2.2, 2.2),
                       (1.2, 1.2, 3.2),
                       (1.2, 2.2, 3.2)], dtype=float)

    triangles = np.array([(0,1,2),
                          (2,1,3),(0,4,2)], dtype=int)
    normals = np.array([[-1, 0, 0],
                        [-1, 0, 0],[-1, 0, 0]], dtype=float)

    return points, triangles, normals

def test_trivial():
    # 16 cube
    levels = 4
    # put triangles onto the 8 cube level
    tri_level = 3
    points, triangles, normals = mk_trivial()
    
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

def test_trivial_deep():
    # 1,024 cube
    levels = 18
    # put triangles onto the 8 cube level
    tri_level = 17
    points, triangles, normals = mk_trivial()
    
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
    
def test_merge_identity():
    # 16 cube
    levels = 4
    # put triangles onto the 8 cube level
    tri_level = 3
    points, triangles, normals = mk_trivial()
    
    tree = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    nulltree = Tree(levels)
    
    s = TriangleSorter.TreeSummer(levels, tri_level)
    s.add(tree)
    s.add(nulltree)
    new = s.finish()
    
    assert trees_with_triIds_equal(new, tree, tri_level)

def test_merge_simple():
    # 16 cube
    levels = 4
    # put triangles onto the 8 cube level
    tri_level = 3
    points, triangles, normals = mk_trivial()
    
    t1 = TriangleSorter.TrianglesToTree_Worker(levels, tri_level, points, triangles[:1], 0)
    t2 = TriangleSorter.TrianglesToTree_Worker(levels, tri_level, points, triangles[1:], 1)
    
    s = TriangleSorter.TreeSummer(levels, tri_level)
    s.add(t1)
    s.add(t2)
    new = s.finish()
    
    tree = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    assert trees_with_triIds_equal(new, tree, tri_level)

def test_merge_self():
    # Actually, t1 + t1 gives an invalid tree - we expect each triangle to be 
    # handled by exactly one process so adding support for ignoring duplicates 
    # is wasteful.
    
    # The final tree should therefore differ from t1
    
    # 16 cube
    levels = 4
    # put triangles onto the 8 cube level
    tri_level = 3
    points, triangles, normals = mk_trivial()
    
    t1 = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    
    s = TriangleSorter.TreeSummer(levels, tri_level)
    s.add(t1)
    s.add(t1)
    tree = s.finish()
    
    # Trees must be structurally identical
    assert tree == t1
    # Trees must have different IDs
    assert not trees_with_triIds_equal(tree, t1, tri_level)
    
def test_parallel():
    levels = 5
    tri_level = 3
    points, triangles, normals = GetSphereNumpy()
    serial = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    
    p1 = TriangleSorter.TrianglesToTreeParallel(levels, tri_level, points, triangles, 1)
    p2 = TriangleSorter.TrianglesToTreeParallel(levels, tri_level, points, triangles, 2)
    p3 = TriangleSorter.TrianglesToTreeParallel(levels, tri_level, points, triangles, 3)
    p4 = TriangleSorter.TrianglesToTreeParallel(levels, tri_level, points, triangles, 4)
    p8 = TriangleSorter.TrianglesToTreeParallel(levels, tri_level, points, triangles, 8)
    assert trees_with_triIds_equal(serial, p1, tri_level)
    assert trees_with_triIds_equal(serial, p2, tri_level)
    assert trees_with_triIds_equal(serial, p3, tri_level)
    assert trees_with_triIds_equal(serial, p4, tri_level)
    assert trees_with_triIds_equal(serial, p8, tri_level)
    
    return
