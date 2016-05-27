#
# Copyright (C) University College London, 2007-2012, all rights reserved.
#
# This file is part of HemeLB and is CONFIDENTIAL. You may not work
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#

import numpy as np
import itertools

def AddIoletsToTree(tree, tri_level, iolets):
    """As for TriangleSorter.TrianglesToTree, but for inlet/outlet plane
    objects."""
    cube_size = 2**tree.levels
    tri_box_size = 2**tri_level
    non_halo_edges = np.arange(0, cube_size, tri_box_size)
    
    # The upper and lower bounds of the regions of influence of the tri_level nodes.
    # +/- 1 because we have to consider the halo.
    lowers = non_halo_edges - 1
    uppers = non_halo_edges + tri_box_size
    
    for io in iolets:
        c = np.array((io.Centre.x, io.Centre.y, io.Centre.z), dtype=float)
        n = np.array((io.Normal.x, io.Normal.y, io.Normal.z), dtype=float)
        r = io.Radius
        
        n2 = n**2 # [nx^2, ny^2, nz^2]
        # 1 - nx^2 = nx^2 + ny^2 + nz^2 - nx^2 = ny^2 + nz^2
        disc_half_size = r*np.sqrt(1.0 - n2)
        
        mins = c - disc_half_size
        maxs = c + disc_half_size
        ranges = map(xrange,
                     uppers.searchsorted(mins), lowers.searchsorted(maxs))

        for box_ind in itertools.product(*ranges):
            node = tree.GetNode(tree.levels - tri_level, np.array(box_ind), relative=True, create=True)
            try:
                node.iolets.append(io)
            except AttributeError:
                node.iolets = [io]
                pass

        