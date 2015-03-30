import numpy as np
import itertools
import multiprocessing
import HemeOct

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
    tree = TrianglesToTree_Worker(n_levels, tri_level, points, triangles, 0)
    
    for node in tree.root.IterDepthFirst(tri_level, tri_level):
        node.triIds = np.array(node.triIds)
        node.triIds.sort()
        continue
    return tree

def TrianglesToTree_Worker(n_levels, tri_level, points, triangles, tri_index_start):
    """The tri_index_start arg is to allow parallelisation of this - it is
    added to the output ID lists. 
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
            triList.append(i_tri + tri_index_start)
            continue
        continue
    return tree

def TrianglesToTreeParallel(n_levels, tri_level, points, triangles, nprocs=None):
    """As for TrianglesToTree, but runs in parallel.
    """
    n_tri = len(triangles)
    if nprocs is None:
        nprocs = multiprocessing.cpu_count()
        
    # tris per proc 
    tpp = float(n_tri) / float(nprocs)
    procs = []
    # answer queue
    q = multiprocessing.Queue()
    
    for i in xrange(nprocs):
        min_ind = int(i*tpp)
        max_ind = min(int((i+1)*tpp), n_tri)
        p = Worker(q, (n_levels, tri_level, points, triangles[min_ind:max_ind], min_ind))
        p.start()
        procs.append(p)
        continue
    
    summer = TreeSummer(n_levels, tri_level)
    # Merge the results together as they come in
    for i in xrange(nprocs):
        summer.add(q.get())
        continue
    
    # Join the child processes
    for p in procs:
        p.join()
        continue
    
    # Actually combine and sort the arrays of IDs
    return summer.finish()

class TreeSummer(object):
    def __init__(self, levels, tri_level):
        self.levels = levels
        self.tri_level = tri_level
        self.tree = HemeOct.Tree(levels)
        return
    
    def add(self, source):
        """Merges the Octrees but only collates the triangle IDs into a list 
        of arrays to avoid lots of redundant copying.
        """
        for src_node in source.root.IterDepthFirst(self.tri_level, self.tri_level):
            new_ids = np.array(src_node.triIds)
            dest_node = self.tree.GetNode(self.tri_level, src_node.offset, create=True)
            try:
                list_of_id_arrays = dest_node.list_of_id_arrays
            except AttributeError:
                list_of_id_arrays = dest_node.list_of_id_arrays = []
                pass
            list_of_id_arrays.append(new_ids)
            
    def finish(self):
        """Finish the merge by concatenating the list of arrays into a single sorted
        one and then deleting this list of arrays.
        """
        for node in self.tree.root.IterDepthFirst(self.tri_level, self.tri_level):
            try:
                list_of_id_arrays = node.list_of_id_arrays
                node.triIds = np.concatenate(list_of_id_arrays)
                node.triIds.sort()
                del node.list_of_id_arrays
            except AttributeError:
                pass
            
        return self.tree
    pass        
        
class Worker(multiprocessing.Process):
    """Simple process subclass to be our parallel worker.
    """
    def __init__(self, ans_queue, args):
        self.ans_queue = ans_queue
        self.args = args
        super(Worker, self).__init__()
        return
    
    def run(self):
        ans = TrianglesToTree_Worker(*self.args)
        self.ans_queue.put(ans)
        return
    pass