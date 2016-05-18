import numpy as np
from .NodeClassifier import NodeClassifier
from .HemeOct import Node, FluidNode, SolidNode
from .Neighbours import neighbours, inverses
from .HitFinder import HitFinder
from .TestResources.sphere import GetSphereNumpy, InsideSphere
from .TriangleSorter import TrianglesToTree
    
from pytest import raises, set_trace

def test_unknown_voxel():
    """If we know nothing about a voxel at all, it should be an error as it
    should never have been created.
    """
    nc = NodeClassifier()
    node = Node(0, np.array((1,1,1)))
    with raises(AttributeError):
        nc(node)

def get_i(*direction):
    return np.all(neighbours == direction, axis=1).nonzero()[0][0]

def test_simple_voxels():
    """Simple solid and fluid cases.
    """
    nc = NodeClassifier()
    
    # Let's assume (1,1,1) is solid and (2,2,2) is fluid
    solid = Node(0, np.array((1,1,1)))
    i = get_i(1,1,1)
    
    solid.intersections = {i: [(0.5, False, 0)]}
    nc(solid)
    assert solid.fluid_count == 0
    
    fluid = Node(0, np.array((2,2,2)))
    fluid.intersections = {inverses[i]: [(0.5, True, 0)]}
    nc(fluid)
    assert fluid.fluid_count == 1

def test_protrusion():
    """Let's consider a sticky-out bit into the fluid.
    Like the link between x and y:
          ^
    x    / \    y   
        /   \   
        0   1    
    """
    nc = NodeClassifier()
    
    i_xy = get_i(1,0,0)
    x = Node(0, np.array((2,2,2)))
    x.intersections = {i_xy: [(0.25, True, 0),
                              (0.75, False, 1)]}    
    nc(x)
    assert x.fluid_count == 1
    
    # And should work with the intersections either way round.
    x = Node(0, np.array((2,2,2)))
    x.intersections = {i_xy: [(0.75, False, 1),
                              (0.25, True, 0)]}
    nc(x)
    assert x.fluid_count == 1
    
    # Now think about y
    y = Node(0, np.array((3,2,2)))
    y.intersections = {inverses[i_xy]: [(0.25, True, 1),
                                        (0.75, False, 0)]}
    nc(y)
    assert y.fluid_count == 1


def test_level_one_planes():
    for z in (1.5, 2.5, 3.5):
        yield level_one_plane, z
        
def level_one_plane(z):
    """Test a node with children.
    
    Let's assume:
        - our children are at (2,2,2)-(3,3,3) 
        - we're near a plane normal to -z (so fluid above, solid below)
    
    For a few value of z (1.5, 2.5, 3.5) to test different cases. 
    """
    parent = Node(1, (2,2,2))
    parent.triIds = np.arange(1)
    
    points = np.array(([ 0, 0, z],
                       [10, 0, z],
                       [ 0,10, z],
                       [10,10, z]), dtype=float)
    triangles = np.array([[0,1,2],
                          [1,2,3]])
    normals = np.array([[0,0,-1]], dtype=float)
    hf = HitFinder(points, triangles, normals)
    
    hf(parent)
    
    nc = NodeClassifier()
    nc(parent)
    
    for vox in parent.IterDepthFirst(0,0):
        vz = vox.offset[2]
        if vz > z:
            assert vox.fluid_count == 1, "Expected voxel to be fluid but wasn't: " + str(vox.offset)
            if vz > z+1:
                assert isinstance(vox, FluidNode)
            else:
                assert isinstance(vox, Node)
                pass
        else:
            assert vox.fluid_count == 0, "Expected voxel to be solid but wasn't: " + str(vox.offset)
            if vz < z-1:
                assert isinstance(vox, SolidNode)
            else:
                assert isinstance(vox, Node)
                pass
            pass
        pass
    
    expected_fluid = int(4*np.ceil(3.0 - z))
    assert parent.fluid_count == expected_fluid, \
        "z: {}; expected: {}' got: {}".format(z, expected_fluid, parent.fluid_count)

def test_sphere():
    levels = 4
    size = 2**levels
    tri_level = 3
    points, triangles, normals = GetSphereNumpy()
    tree = TrianglesToTree(levels, tri_level, points, triangles)
    
    finder = HitFinder(points, triangles, normals)
    for tn in tree.IterDepthFirst(tri_level, tri_level):
        finder(tn)
    
    nc = NodeClassifier()
    nc(tree.root)
    inside = InsideSphere(np.mgrid[:size,:size,:size].transpose((1,2,3,0)))
    
    # Check the leaf voxels
    for vox in tree.IterDepthFirst(0,0):
        try:
            if inside[tuple(vox.offset)]:
                assert vox.fluid_count == 1, "Voxel should be fluid"
            else:
                assert vox.fluid_count == 0, "Voxel should be solid"
        except AssertionError:
            set_trace()