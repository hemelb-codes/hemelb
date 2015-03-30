from .HemeOct import Tree
from .Oct import NodeError
import numpy as np

def mk_sphere_tree():
    pwr = 4
    c = np.array((7.5, 7.5, 7.5))
    r = 6.0
    def dist(xyz):
        return np.sqrt(np.sum((xyz - c)**2, axis=-1))
    
    tree = Tree(pwr)
    size = 2**pwr
    for ijk in np.ndindex((size, size, size)):
        if dist(ijk) <= r:
            node = tree.GetNode(0, np.array(ijk), create=True)
            node.fluid_count = 1
    return tree

class test_withSphereTree(object):
    """Make a 16 cube tree with a sphere radius = 6 of fluid centred on
    (7.5, 7.5, 7.5).
    
    All fluid leaf nodes exist.
    """
    def setup(self):
        self.tree = mk_sphere_tree()
        return
    
    def test_init(self):
        i,j,k = np.mgrid[:16,:16,:16]
        x = i - 7.5
        y = j - 7.5
        z = k - 7.5
        r = np.sqrt(x**2 + y**2 + z**2)
        
        for ijk, rval in np.ndenumerate(r):
            if rval <= 6:
                n = self.tree.GetNode(0, np.array(ijk))
                assert n.fluid_count == 1
            else:
                try:
                    n = self.tree.GetNode(0, np.array(ijk))
                except NodeError:
                    pass
                else:
                    raise ValueError("Node existed that shouldn't: " + ijk)
                pass
            continue
        return
    