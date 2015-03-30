import numpy as np
import itertools
import multiprocessing
import HemeOct
import pdb

def TrianglesToTree(n_levels, tri_level, points, triangles):
    """Create an Octree with n_levels potential number of levels.
    Nodes will be created down to the level tri_level but only if voxel-voxel
    links originating within that node could potentially intersect the triangles
    specified by (points, triangles), where:
    
    points is a numpy array of doubles with shape == (nPoints, 3) specifying
        the vertices.
        
    triangles is a numpy array of ints with shape == (nTris, 3) specifiying the
        IDs of the verticies comprising the triangle.
        
    The create nodes at tri_level will have an attribute 'triIds' added which is
    a list of the IDs of the triangles that might intersect its voxel-voxel links.
    """
    cube_size = 2**n_levels
    tri_box_size = 2**tri_level
    non_halo_edges = np.arange(0, cube_size, tri_box_size)
    
    # The upper and lower bounds of the regions of influence of the tri_level nodes.
    # +/- 1 because we have to consider the halo.
    lowers = non_halo_edges - 1
    uppers = non_halo_edges + tri_box_size
    
    tree = HemeOct.Tree(n_levels)
    for i_tri, triPtIds in enumerate(triangles):
        # bounding box of the triangle
        mins = points[triPtIds].min(axis=0)
        maxs = points[triPtIds].max(axis=0)
        ranges = map(xrange,
                     uppers.searchsorted(mins), lowers.searchsorted(maxs))
        # This bit is the only non-thread-safe part.
        # We will later merge the trees if this was run in parallel
        for box_ind in itertools.product(*ranges):
            node = tree.GetNode(n_levels - tri_level, np.array(box_ind), relative=True, create=True)
            try:
                triList = node.triIds
            except AttributeError:
                triList = node.triIds = []
                pass
            triList.append(i_tri)
            continue
        continue
    
    for node in tree.root.IterDepthFirst(tri_level, tri_level):
        node.triIds = np.array(node.triIds)
        node.triIds.sort()
        continue
    return tree
