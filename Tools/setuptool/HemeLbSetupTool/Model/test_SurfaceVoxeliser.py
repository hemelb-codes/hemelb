import numpy as np
from TestResouces.simple_meshes import mk_trivial
from TestResources.sphere import GetSphereNumpy, Radius
from .SurfaceVoxeliser import SurfaceVoxeliser
from .HemeOct import Tree
import Oct

from pytest import set_trace
def trees_with_triIds_equal(t1, t2, tri_level):
    if t1 != t2:
        return False
    # If here, they're structurally the same, check triIds
    for n1 in t1.IterDepthFirst(tri_level, tri_level):
        n2 = t2.GetNode(n1.levels, n1.offset)
        if not n1.triIds == n2.triIds:
            return False
        continue
    # None of the nodes differ so we are true!
    return True

def test_trivial_points():
    # 8 cube
    levels = 3
    n = 2**levels
    points, triangles, normals = mk_trivial()
    voxer = SurfaceVoxeliser(points, triangles, normals, levels)
    
    # Get full grid of points
    coords = np.mgrid[:n,:n,:n].reshape((3,n**3)).transpose()
    mask = np.empty(n**3, dtype=bool)
    
    for iPt in xrange(len(points)):
        mask[:] = False
        voxer.FilterPoint(iPt, coords, mask)
    
        # innies
        dist2 = np.sum((coords[mask] - points[iPt])**2,axis=1)
        assert np.all(dist2 <= 3.0/4.0)
    
        # outies
        dist2 = np.sum((coords[np.logical_not(mask)] - points[iPt])**2,axis=1)
        assert np.all(dist2 > 3.0/4.0)

def test_trivial_edges():
    levels = 3
    n = 2**levels
    points, triangles, normals = mk_trivial()
    voxer = SurfaceVoxeliser(points, triangles, normals, levels)
    
    # Get full grid of points
    coords = np.mgrid[:n,:n,:n].reshape((3,n**3)).transpose()
    mask = np.empty(n**3, dtype=bool)
    
    for n, (i,j) in enumerate(((0,1), (1,2), (2,0))):
        mask[:] = False
        voxer.FilterEdge(i, j, coords, mask)
        if n==0:
            # cylinder along z axis
            
            # innies
            x,y,z = coords[mask].transpose()
            assert np.all(z >= 1.2)
            assert np.all(z <= 2.2)
            r2 = (x-1.2)**2 + (y-1.2)**2
            assert np.all(r2 <= 3.0/4.0)
            
            # outies
            x,y,z = coords[np.logical_not(mask)].transpose()
            out = z < 1.2
            out |= z > 2.2
            r2 = (x-1.2)**2 + (y-1.2)**2
            out |= r2 > 3.0/4.0
            assert np.all(out)
            
        if n == 2:
            # y cylinder
            # innies
            x,y,z = coords[mask].transpose()
            assert np.all(y >= 1.2)
            assert np.all(y <= 2.2)
            r2 = (x-1.2)**2 + (z-1.2)**2
            assert np.all(r2 <= 3.0/4.0)
            
            # outies
            x,y,z = coords[np.logical_not(mask)].transpose()
            out = y < 1.2
            out |= y > 2.2
            r2 = (x-1.2)**2 + (z-1.2)**2
            out |= r2 > 3.0/4.0
            assert np.all(out)
            

def test_trivial_plane():
    levels = 3
    n = 2**levels
    points, triangles, normals = mk_trivial()
    voxer = SurfaceVoxeliser(points, triangles, normals, levels)
    
    # Get full grid of points
    coords = np.mgrid[:n,:n,:n].reshape((3,n**3)).transpose()
    mask = np.zeros(n**3, dtype=bool)
    
    # This tri shouldn't get any points
    voxer.FilterTriangle(0, coords, mask)
    assert np.all(np.logical_not(mask))
    
    # This one should get (1,2,2)
    voxer.FilterTriangle(1, coords, mask)
    pt122_mask = np.all(coords == (1,2,2), axis=1)
    
    assert np.all(mask[pt122_mask])
    assert np.all(np.logical_not(mask[np.logical_not(pt122_mask)]))
    
def test_trivial():
    # 16 cube
    levels = 4
    points, triangles, normals = mk_trivial()

    voxer = SurfaceVoxeliser(points, triangles, normals, levels)
    voxer.Execute()
    tree = voxer.Tree    
    
    assert tree.levels == levels
    
    expected_nodes = set(((1,1,1),
                          (1,1,2),
                          (1,2,1),
                          (1,2,2),
                          
                          (2,1,1),
                          (2,1,2),
                          (2,2,1),
                          (2,2,2),
                          
                          (1,1,3),
                          (1,2,3),
                          (1,3,1),
                          (1,3,2)))
    
    for node in tree.IterDepthFirst():
        if node.levels == 0:
            # leaf
            ijk = tuple(node.offset)
            # this will raise KeyError if not present
            expected_nodes.remove(ijk)
        elif node.levels == 1:
            assert np.all(np.logical_or(node.offset == 0, node.offset == 2))
        else:
            # All higher levels only in the zero octant
            assert np.all(node.offset == 0)
            
        continue
    

def overlap1d(amin, amax, bmin, bmax):
    return np.logical_and(amin < bmax, amax > bmin)
def overlap3d(amin, amax, bmin, bmax):
    return np.logical_and(
        np.logical_and(
            overlap1d(amin[..., 0], amax[..., 0], bmin[..., 0], bmax[..., 0]),
            overlap1d(amin[..., 1], amax[..., 1], bmin[..., 1], bmax[..., 1])),
        overlap1d(amin[..., 2], amax[..., 2], bmin[..., 2], bmax[..., 2]))
    
def test_sphere_tri90():
    # Initial testing showed that this fails for some triangles, including this
    levels = 5
    n = 2**levels
    points, triangles, normals = GetSphereNumpy()
    voxer = SurfaceVoxeliser(points, triangles, normals, levels)
    
    iTri = 90
    tri_pt_ids = triangles[iTri]
    tri_pts = points[tri_pt_ids]
    norm = normals[iTri]
    
    lo, hi = voxer.AABB(iTri)
    lo = lo.astype(int) - 1
    # 1 to round up, 1 to be safe, 1 to get range upper bound
    hi = hi.astype(int) + 3
    
    # This region of interest bounds the triangle safely
    roi = np.mgrid[lo[0]:hi[0],
                   lo[1]:hi[1],
                   lo[2]:hi[2]]
    
    voxels = roi.reshape((3, np.prod(hi-lo))).transpose()
    inside_mask = np.zeros(len(voxels), dtype=bool)
    union = np.zeros(len(voxels), dtype=bool)
    
    expected_voxels = []
    
    ev = set()
    ev.add((22,8,18))
    ev.add((22,9,18))
    ev.add((23,9,18))
    expected_voxels.append(ev)
    
    ev = set()
    ev.add((22,8,13))
    ev.add((22,9,13))
    ev.add((23,9,13))
    expected_voxels.append(ev)
    
    ev = set()
    ev.add((25,15,13))
    ev.add((25,16,13))
    expected_voxels.append(ev)

    for i, iPt in enumerate(tri_pt_ids):
        inside_mask[:] = False
        voxer.FilterPoint(iPt, voxels, inside_mask)
        ev = expected_voxels[i]
        for v in voxels[inside_mask]:
            v = tuple(v.tolist())
            # raises KeyError if no there
            ev.remove(v)
        assert len(ev) == 0
        union |= inside_mask
    
    # Can't face working this out by hand!
    for i in xrange(3):
        inside_mask[:] = False
        iPt = tri_pt_ids[i]
        jPt = tri_pt_ids[(i+1)%3]
        voxer.FilterEdge(iPt, jPt, voxels, inside_mask)
        union |= inside_mask
    
    # or this...
    inside_mask[:] = False
    voxer.FilterTriangle(iTri, voxels, inside_mask)
    union |= inside_mask
        

def connected_region(array, idx):
    """Helper function that does a 6-connected flood fill to mark
    the continuous region with the same value as idx.
    """
    nx, ny, nz = array.shape
    
    ans = np.zeros_like(array, dtype=bool)        
    checked = np.zeros_like(ans)
    target = array[idx]
    
    dx = (1, 0, 0,-1, 0, 0)
    dy = (0, 1, 0, 0,-1, 0)
    dz = (0, 0, 1, 0, 0,-1)
    
    stack = []
    
    stack.append(idx)
    while len(stack):
        x,y,z = stack.pop()
        ans[x,y,z] = True
        checked[x,y,z] = True
        for i in xrange(6):
            u = x+dx[i]
            if u < 0 or u >= nx:
                continue
            v = y+dy[i]
            if v < 0 or v >= ny:
                continue
            w = z+dz[i]
            if w < 0 or w >= nz:
                continue
            
            if not checked[u,v,w] and array[u,v,w] == target:
                stack.append((u,v,w))
                pass
            continue
        continue
    
    return ans
def test_connected1():
    im = np.array((
                   ((0,0,0,0,0),
                    (0,1,1,1,0),
                    (0,1,1,1,0),
                    (0,1,1,1,0),
                    (0,0,0,0,0),
                   ),
                   ((0,0,0,0,0),
                    (0,1,1,1,0),
                    (0,1,1,1,0),
                    (0,1,1,1,0),
                    (0,0,0,0,0),
                    ),
                   ((0,0,0,0,0),
                    (0,1,1,1,0),
                    (0,1,1,1,0),
                    (0,1,1,1,0),
                    (0,0,0,0,0),
                    ),
                   ), dtype=bool)
    ff = connected_region(im, (0,2,2))
    assert np.all(ff == im)
    
def test_connected2():
    im = np.array((
                   ((1,1,0,0,0),
                    (1,1,0,0,0),
                    (0,0,0,0,0),
                    (0,0,0,1,1),
                    (0,0,0,1,1),
                   ),
                   ), dtype=bool)
    ff = connected_region(im, (0,0,0))
    assert np.all(ff == np.array((
                                  ((1,1,0,0,0),
                                   (1,1,0,0,0),
                                   (0,0,0,0,0),
                                   (0,0,0,0,0),
                                   (0,0,0,0,0),
                                   ),
                                  ), dtype=bool))

def test_sphere():
    levels = 5
    n = 2**levels
    points, triangles, normals = GetSphereNumpy()
    voxer = SurfaceVoxeliser(points, triangles, normals, levels)
    voxer.Execute()
    tree = voxer.Tree

    edge_mask = Oct.TreeToMaskArray(tree)
    
    interior = connected_region(edge_mask, (15,15,15))
    x,y,z = interior.nonzero()
    x = x - 15.5
    y = y - 15.5
    z = z - 15.5
    r2 = x**2 + y**2 + z**2
    assert np.all(r2 < 100.0)
    
    seen_tris = set()
    
    for node in tree.IterDepthFirst(0,0):
        assert Radius(node.offset) <= 10.75
        seen_tris.update(node.triIds)
        continue
    
    seen_tris = list(seen_tris)
    seen_tris.sort()
    nTri = len(triangles)
    assert np.all(np.arange(nTri) == seen_tris)
    
def test_tovtk():
    levels = 5
    points, triangles, normals = GetSphereNumpy()
    voxer = SurfaceVoxeliser(points, triangles, normals, levels)
    voxer.Execute()
    tree = voxer.Tree
    tree.ToVtk('sphere.vtk')
