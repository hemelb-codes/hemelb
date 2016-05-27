#
# Copyright (C) University College London, 2007-2012, all rights reserved.
#
# This file is part of HemeLB and is CONFIDENTIAL. You may not work
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#
import numpy as np
from TestResources.simple_meshes import mk_trivial, mk_trivial2
from .IoletSorter import AddIoletsToTree
from .HemeOct import Tree
from .Iolets import Inlet, Outlet
from .Vector import Vector

from pytest import set_trace
def mk_inlet_z():
    inlet = Inlet(Centre=Vector(4.5, 2.4, 2.8),
                  Normal=Vector(0.0, 0.0, 1.0),
                  Radius=2.0)
    return inlet

def mk_inlet_squint():
    n = np.array((2, -7, 11), dtype=float)
    n /= np.sum(n**2)
    return Inlet(Centre=Vector(8,9,10),
                 Normal=Vector(*n),
                 Radius=3.9)
    
def test_simple():
    # 16 cube
    levels = 4
    tri_level = 2
    tree = Tree(levels)
    iolets = [mk_inlet_z()]
    AddIoletsToTree(tree, tri_level, iolets)
    
    offsets = set(((0, 0, 0),
                   (0, 4, 0),
                   (4, 0, 0),
                   (4, 4, 0)))
    for node in tree.IterDepthFirst(tri_level, tri_level):
        offsets.remove(tuple(node.offset.tolist()))
    
    assert(len(offsets) == 0)
    
def test_squint():
    # 16 cube
    levels = 4
    tri_level = 2
    tree = Tree(levels)
    iolets = [mk_inlet_squint()]
    AddIoletsToTree(tree, tri_level, iolets)
    
    #8,9,10 - 3.9
    offsets = set(((4,4,4),
                   (4,4,8),
                   (4,4,12),
                   (4,8,4),
                   (4,8,8),
                   (4,8,12),
                   (4,12,4),
                   (4,12,8),
                   (4,12,12),
                   (8,4,4),
                   (8,4,8),
                   (8,4,12),
                   (8,8,4),
                   (8,8,8),
                   (8,8,12),
                   (8,12,4),
                   (8,12,8),
                   (8,12,12),
                   (12,4,4),
                   (12,4,8),
                   (12,4,12),
                   (12,8,4),
                   (12,8,8),
                   (12,8,12),
                   (12,12,4),
                   (12,12,8),
                   (12,12,12)
                   ))
                    
    for node in tree.IterDepthFirst(tri_level, tri_level):
        offsets.remove(tuple(node.offset.tolist()))
    
    assert(len(offsets) == 0)
    