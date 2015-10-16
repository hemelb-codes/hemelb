import numpy as np
from .HitFinder import HitFinder, Face, FaceGroup
from .TriangleSorter import TrianglesToTree
from .test_TriangleSorter import mk_trivial, mk_trivial2
from .Neighbours import neighbours, inverses
from TestResources.sphere import GetSphereNumpy, InsideSphere
from . import HemeOct
from . import Oct

def approx_eq(x, y, tol):
    return np.all(np.abs(np.subtract(x,y)) < tol)
def approx_lt(x, y, tol):
    return x < (y + tol)

def test_face():
    xy = Face(4, (0,0,1))
    assert xy.name == '+xy'
    assert xy.iter_dims == [0,1]
    seen = np.zeros((6,6), dtype=bool)
    for ijk in xy:
        assert ijk[2] == -1
        seen[ijk[0], ijk[1]] = True
        continue
    assert np.all(seen)
    
    zy = Face(5, (-1, 0, 0))
    assert zy.name == '-yz'
    assert zy.iter_dims == [1,2]
    
    seen = np.zeros((7,7), dtype=bool)
    for ijk in zy:
        assert ijk[0] == 5
        seen[tuple(ijk[1:]+1)] = True
        continue
    assert np.all(seen)
    
    xz = Face(5, (0,1,0))
    assert xz.name == '+xz'
    assert xz.iter_dims == [0,2]
    
    seen = np.ones((7,7,7), dtype=bool)
    seen[6,:,:] = 0
    seen[:,0,:] = 0
    for ijk in zy:
        seen[tuple(ijk+1)] ^= 1
    for ijk in xz.Iter(zy):
        seen[tuple(ijk+1)] ^= 1
    assert np.all(seen)
    
    xy = Face(5, (0,0,1))
    assert xy.name == '+xy'
    assert xy.iter_dims == [0,1]
    
    seen = np.ones((7,7,7), dtype=bool)
    seen[6,:,:] = 0
    seen[:,0,:] = 0
    seen[:,:,0] = 0
    for ijk in zy:
        seen[tuple(ijk+1)] ^= 1
    for ijk in xz.Iter(zy):
        seen[tuple(ijk+1)] ^= 1
    for ijk in xy.Iter(zy, xz):
        seen[tuple(ijk+1)] ^= 1
    assert np.all(seen)
    
def test_face_group():
    zy = Face(5, (-1, 0, 0))
    xz = Face(5, (0,1,0))
    xy = Face(5, (0,0,1))
    grp = FaceGroup(xy, xz, zy)
    seen = np.ones((7,7,7), dtype=bool)
    seen[6,:,:] = 0
    seen[:,0,:] = 0
    seen[:,:,0] = 0
    for ijk in grp:
        seen[tuple(ijk+1)] ^= 1
    assert np.all(seen)

def test_ray_list():
    # Monkey patch the finder with a DoRay method that just records what it was
    # called with
    finder = HitFinder([], [], [])
    def _DoRay(tri_node, i_vec, startIjk, endIjk):
        vec = neighbours[i_vec]
        i_opp = inverses[i_vec]
        ijk = startIjk + vec
        while not np.all(ijk == endIjk):
            seen[i_vec][tuple(ijk - tri_node.offset)] ^= 1
            seen[i_opp][tuple(ijk - tri_node.offset)] ^= 1
            ijk += vec
            continue
    finder._DoRay = _DoRay
    
    tri_node = HemeOct.Node(1, (0,0,0))
    seen = np.zeros((26, 2,2,2), dtype=bool)
    finder(tri_node)
    assert np.all(seen)
    
    tri_node = HemeOct.Node(2, (4, 0, 8))
    seen = np.zeros((26, 4,4,4), dtype=bool)
    finder(tri_node)
    assert np.all(seen)
    
def test_doray1():
    points, triangles, normals = mk_trivial()
    node = HemeOct.Node(2, (0,0,0))
    i_x = 21
    assert np.all(neighbours[i_x] == (1,0,0))
    
    finder = HitFinder(points, triangles, normals)    
    for y,z in np.ndindex(4,4):
        finder._DoRay(node, i_x, np.array((-1,y,z,)), np.array((4,y,z)))
        if y == 2 and z == 2:
            n1  = node.GetNode(0, np.array((1,y,z)))
            n2  = node.GetNode(0, np.array((2,y,z)))
            # 1 direction with hits each
            assert len(n1.intersections) == 1
            assert len(n2.intersections) == 1
            
            hits1 = n1.intersections[i_x]
            hits2 = n2.intersections[inverses[i_x]]
            # 1 hit per link
            assert len(hits1) == 1
            assert len(hits2) == 1
            hit1 = hits1[0]
            hit2 = hits2[0]
            # the t's total 1
            assert np.abs(hit1[0] + hit2[0] - 1.0) < 1e-9
            # opposite directions
            assert hit1[1] != hit2[1]
            # same tri
            assert hit1[2] == hit2[2]
            
        else:
            for x in xrange(4):
                missing = False
                try:
                    node.GetNode(0, np.array((x,y,z)))
                except Oct.NodeError:
                    missing = True
                    pass
                assert missing, "Unexpected node: "+str((x,y,z))
                continue
            pass
        continue

def test_doray2():
    # As above, but there should be no hits at all
    points, triangles, normals = mk_trivial()
    node = HemeOct.Node(2, (0,0,0))
    i_z = 13
    assert np.all(neighbours[i_z] == (0,0,1))
    
    finder = HitFinder(points, triangles, normals)
    for x,y in np.ndindex(4,4):
        finder._DoRay(node, i_z, np.array((x,y,-1)), np.array((x,y,4)))
        for z in xrange(4):
            missing = False
            try:
                node.GetNode(0, np.array((x,y,z)))
            except Oct.NodeError:
                missing = True
                assert missing, "Unexpected node: "+str((x,y,z))
            continue
        continue
    
def test_square():
    points, triangles, normals = mk_trivial()
    tree = HemeOct.Tree(2)
    root = tree.root
    finder = HitFinder(points, triangles, normals)
    
    finder(root)
    nhits = 0
    for vox in root.IterDepthFirst(0,0):
        
        for i_vec, hits in vox.intersections.iteritems():
            vec = neighbours[i_vec]
            neigh_offset = vox.offset + vec
            neigh = root.GetNode(0, neigh_offset)
            opp = inverses[i_vec]
            # There must be a hit coming in the other direction
            assert opp in neigh.intersections, "Neighbour must have a hit for the reverse link"
            opp_hits = neigh.intersections[opp]
            assert len(hits) == len(opp_hits), "Must have the same number of hits"
            
            # Sort them based on the t, the distance along the link
            hits.sort(key=lambda x: x[0])
            opp_hits.sort(key=lambda x: -x[0])
            n = len(hits)
            nhits += n
            for i in xrange(n):
                assert approx_eq(hits[i][0] + opp_hits[i][0], 1.0, 1e-8), "The lengths along the vector must sum to 1"
                assert hits[i][1] != opp_hits[i][1], "The normal flags must be opposite"
                assert hits[i][2] == opp_hits[i][2], "Must hit the same triangle"
                t = hits[i][0]
                hit_point = vox.offset + t * vec
                
                assert approx_eq(hit_point[0], 1.2, 1e-9), "Hit point coord wrong"
                assert approx_lt(1.2, hit_point[1], 1e-9) and approx_lt(hit_point[1], 2.2, 1e-9)
                assert approx_lt(1.2, hit_point[2], 1e-9) and approx_lt(hit_point[2], 2.2, 1e-9)
                
        continue
    assert nhits == 32


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
            opp_hits.sort(key=lambda x: -x[0])
            n = len(hits)
            for i in xrange(n):
                assert approx_eq(hits[i][0] + opp_hits[i][0], 1.0, 1e-8), "The lengths along the vector must sum to 1"
                assert hits[i][1] != opp_hits[i][1], "The normal flags must be opposite"
                assert hits[i][2] == opp_hits[i][2], "Must hit the same triangle"
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
