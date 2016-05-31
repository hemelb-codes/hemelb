import numpy as np
import TriangleSorter
from SurfaceVoxeliser import SurfaceVoxeliser
from VoxelClassifier import VoxelClassifier
from . import Oct
from TestResources.simple_meshes import mk_trivial
from TestResources.sphere import GetSphereNumpy, Radius
from pytest import set_trace

def test_trivial():
    # 16 cube
    levels = 4
    tri_level = 2
    points, triangles, normals = mk_trivial()
    labels = -np.ones(len(triangles), dtype=int)
    
    tree = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    voxer = SurfaceVoxeliser(points, triangles, normals, tree, tri_level)
    voxer.Execute()
    tree = voxer.Tree
    
    vc = VoxelClassifier(points, triangles, normals, labels)
    new_tree = vc(tree, tri_level)
    
    for voxel in new_tree.IterDepthFirst(0,0):
        # Only nodes with x = 2 should exist
        assert voxel.offset[0] == 2
        # The cut distances should be ~0.8
        cut_mask = voxel.cut_tri >= 0
        real_cd = voxel.cut_dist[cut_mask]
        assert np.all((real_cd - 0.8)**2 < 1e-6)
    return

def test_sphere():
    levels = 5
    tri_level = 3

    points, triangles, normals = GetSphereNumpy()
    labels = -np.ones(len(triangles), dtype=int)
    
    tree = TriangleSorter.TrianglesToTree(levels, tri_level, points, triangles)
    voxer = SurfaceVoxeliser(points, triangles, normals, tree, tri_level)
    voxer.Execute()
    tree = voxer.Tree
    
    vc = VoxelClassifier(points, triangles, normals, labels)
    new_tree = vc(tree, tri_level)
    for voxel in new_tree.IterDepthFirst(0,0):
        r = Radius(voxel.offset)
        assert r < 10.0
        
    fluid_edge_mask = Oct.TreeToMaskArray(new_tree)
    np.save("fluid_edge_sphere.npy", fluid_edge_mask)