from .Oct import Tree, NodeError
import numpy as np
from pytest import raises

def test_Octree_create():
    tree = Tree(3)
    # Get the root relatively
    root1 = tree.GetNode(0, np.array((0,0,0)), relative=True)
    # Now absolutely
    root2 = tree.GetNode(3, np.array((0,0,0)))
    assert root1 is root2
    assert root1.children is None
    assert tree.levels == 3

def mk_8cube_234():
    # represent an 8x8x8 cube with only the leaf at 2,3,4 existing
    tree = Tree(3)
    tree.GetNode(0, np.array((2,3,4)), create=True)
    return tree

def test_Octree_simple():
    def check_child(node, ind):
        for ijk, child in np.ndenumerate(node.children):
            if ijk == tuple(ind):
                assert child is not None
            else:
                assert child is None

    tree = mk_8cube_234()
    # (2,3,4) in binary is:
    # (010, 011, 100)
    # which tells us the path down the tree to the leaf
    n0 = tree.root
    assert n0.levels == 3
    
    ind = np.array((0,0,1))
    check_child(n0, ind)
    
    n1 = tree.GetNode(1, ind, relative=True)
    assert np.all(n1.rel_offset == ind) 
    assert np.all(n1.offset == ind*4)
    n1_abs = tree.GetNode(2, ind*2**2)
    assert n1 is n1_abs
    
    ind = np.array((1,1,0))
    check_child(n1, ind)
    
    n2 = n1.GetNode(1, ind, relative=True)
    assert np.all(n2.offset == (2,2,4))
    assert np.all(n2.rel_offset == (1,1,2))
    
    n2_rel = tree.GetNode(2, np.array((1,1,2)), relative=True)
    n2_abs = tree.GetNode(1, np.array((2,2,4)))
    assert n2 is n2_rel
    assert n2 is n2_abs
    
    ind = np.array((0,1,0))
    check_child(n2, ind)
    
    n3 = n2.GetNode(1, ind, relative=True)
    assert np.all(n3.offset == (2,3,4))
    assert np.all(n3.rel_offset == (2,3,4))
    
    n3_rel = tree.GetNode(3, np.array((2,3,4)), relative=True)
    n3_abs = tree.GetNode(0, np.array((2,3,4)))
    assert n3 is n3_rel
    assert n3 is n3_abs
    
    assert n3.children is None
    
def test_get_nonexistant():
    tree = mk_8cube_234()
    with raises(NodeError):
        tree.GetNode(0, np.array((0,0,0)))
    
def test_iter():
    tree = mk_8cube_234()
    expected_nodes = np.array([[2,3,4],
                               [2,2,4],
                               [0,0,4],
                               [0,0,0]])
    # Do the whole tree
    i = 0
    for node in tree.IterDepthFirst():
        assert node.levels == i
        assert np.all(node.offset == expected_nodes[i])
        i += 1    
    assert i == 4
    
    # Go down to the penultimate level
    i = 1
    for node in tree.IterDepthFirst(1):
        assert node.levels == i
        assert np.all(node.offset == expected_nodes[i])
        i += 1    
    assert i == 4
    
    # Only 1 & 2
    i = 1
    for node in tree.IterDepthFirst(1,2):
        assert node.levels == i
        assert np.all(node.offset == expected_nodes[i])
        i += 1
    assert i == 3
    
def test_delnode():
    tree = mk_8cube_234()
    
    tree.DelNode(0, np.array((2,3,4)))
    expected_nodes = np.array([[2,3,4],
                               [2,2,4],
                               [0,0,4],
                               [0,0,0]])
    # Do the whole tree
    i = 1 # Cos we deleted the leaf
    for node in tree.IterDepthFirst():
        assert node.levels == i
        assert np.all(node.offset == expected_nodes[i])
        i += 1    
    assert i == 4
    return

def test_eq():
    t1 = mk_8cube_234()
    t2 = mk_8cube_234()
    
    assert t1 == t2
    t3 = Tree(3)
    assert t1 != t3
    
    t2.DelNode(0, np.array((2,3,4)))
    assert t1 != t2
    t2.DelNode(1, np.array((0,0,1)), relative=True)
    assert t2 == t3
    