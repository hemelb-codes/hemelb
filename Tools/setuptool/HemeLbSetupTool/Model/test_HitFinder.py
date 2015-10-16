import numpy as np
from .HitFinder import HitFinder
from .TriangleSorter import TrianglesToTree
from .test_TriangleSorter import mk_trivial, mk_trivial2
from .Neighbours import neighbours, inverses
from TestResources.sphere import GetSphereNumpy, InsideSphere

def approx_eq(x, y, tol):
    return np.all(np.abs(np.subtract(x,y)) < tol)

def test_square():
    # 16 cube
    levels = 4
    # put triangles onto the 8 cube level
    tri_level = 3

    points, triangles, normals = mk_trivial()
    tree = TrianglesToTree(levels, tri_level, points, triangles)
    
    tri_node = tree.GetNode(tri_level, np.array([0,0,0]))
    finder = HitFinder(points, triangles, normals)
    
    finder(tri_node)
    n = 0
    for vox in tri_node.IterDepthFirst(0,0):
        
        for i_vec, hits in vox.intersections.iteritems():
            vec = neighbours[i_vec]
            neigh_offset = vox.offset + vec
            neigh = tri_node.GetNode(0, neigh_offset)
            opp = inverses[i_vec]
            # There must be a hit coming in the other direction
            assert opp in neigh.intersections, "Neighbour must have a hit for the reverse link"
            opp_hits = neigh.intersections[opp]
            assert len(hits) == len(opp_hits), "Must have the same number of hits"
            
            # Sort them based on the t, the distance along the link
            hits.sort(key=lambda x: x[0])
            opp_hits.sort(key=lambda x: -x[0])
            n = len(hits)
            for i in xrange(n):
                assert approx_eq(hits[i][0] + opp_hits[i][0], 1.0, 1e-8), "The lengths along the vector must sum to 1"
                assert hits[i][1] != opp_hits[i][1], "The normal flags must be opposite"
                assert hits[i][2] == opp_hits[i][2], "Must hit the same triangle"
            
        n += 1
        continue
    assert n == 13


def test_rectangle():
    # 16 cube
    levels = 4
    # put triangles onto the 8 cube level
    tri_level = 3

    points, triangles, normals = mk_trivial2()
    tree = TrianglesToTree(levels, tri_level, points, triangles)
    
    tri_node = tree.GetNode(tri_level, np.array([0,0,0]))
    finder = HitFinder(points, triangles, normals)
    
    finder(tri_node)
    m = 0
    for vox in tri_node.IterDepthFirst(0,0):
        
        for i_vec, hits in vox.intersections.iteritems():
            vec = neighbours[i_vec]
            neigh_offset = vox.offset + vec
            neigh = tri_node.GetNode(0, neigh_offset)
            opp = inverses[i_vec]
            # There must be a hit coming in the other direction
            assert opp in neigh.intersections, "Neighbour must have a hit for the reverse link"
            opp_hits = neigh.intersections[opp]
            assert len(hits) == len(opp_hits), "Must have the same number of hits"
            
            # Sort them based on the t, the distance along the link
            hits.sort(key=lambda x: x[0])
            opp_hits.sort(key=lambda x: x[0])
            n = len(hits)
            for i in xrange(n):
                assert approx_eq(hits[i][0] + opp_hits[n-1-i][0], 1.0, 1e-8), "The lengths along the vector must sum to 1"
                assert hits[i][1] != opp_hits[n-1-i][1], "The normal flags must be opposite"
                assert hits[i][2] == opp_hits[n-1-i][2], "Must hit the same triangle"
        m += 1
        continue
    
def test_sphere():
    levels = 5
    size = 2**levels
    tri_level = 3
    points, triangles, normals = GetSphereNumpy()
    tree = TrianglesToTree(levels, tri_level, points, triangles)
    
    finder = HitFinder(points, triangles, normals)
    for tn in tree.IterDepthFirst(tri_level, tri_level):
        finder(tn)
    
    inside = InsideSphere(np.mgrid[:size,:size,:size].transpose((1,2,3,0)))
    
    for ijk, fluid in np.ndenumerate(inside):
        for i_vec, vec in enumerate(neighbours):
            neigh_ijk = ijk+vec
            if np.any(neigh_ijk < 0) or np.any(neigh_ijk >= size):
                continue
            if inside[tuple(neigh_ijk)] != fluid:
                # There ought to be a hit here!
                try:
                    node = tree.GetNode(0, np.array(ijk))
                    hits = node.intersections[i_vec]
                    hits.sort(key=lambda h: h[0])
                    assert hits[0][1] == fluid
                except:
                    from nose.tools import set_trace
                    set_trace()
