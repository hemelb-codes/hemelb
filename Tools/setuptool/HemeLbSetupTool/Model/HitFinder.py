import numpy as np
from Neighbours import neighbours

import pdb
class HitFinder(object):
    """This class's job is to take a triangle node (from TriangleSorter) and 
    add all the intersections between the triangles and its contained voxels. 
    
    Each voxel node with intersections will be created and get an attribute
    'intersections' which is a map from the vector ID (see Neighbours) to a list
    of the hits. A hit is a tuple of (distance_along_vector, outward_flag, tri_ID).
    
    """
    def __init__(self, points, triangles, normals):
        self.triangles = triangles
        self.points = points
        self.normals = normals
        self.a = points[triangles[:,0]]
        self.b = points[triangles[:,1]]
        self.ab = self.b - self.a
        self.c = points[triangles[:,2]]
        self.ac = self.c - self.a

        return
        
    def __call__(self, tri_node):
        triIds = tri_node.triIds
        del tri_node.triIds
        aLocal = self.a[triIds]
        nLocal = self.normals[triIds]
        aDOTn = np.sum(aLocal * nLocal, 1)

        size = 2**tri_node.levels

        pGrid = np.mgrid[:size,:size,:size]
        pGrid = pGrid.reshape((3, size**3)).transpose()
        pGrid += tri_node.offset
        pDOTn = np.tensordot(pGrid, nLocal, [[1], [1]])
        
        total_hits = 0
        
        for iDisp, disp in enumerate(neighbours):
            # Need to figure out if the line from pGrid[i] to qGrid[i]  = p + disp intersects the surface
            # First, does it intersect the plane defined by a and the normal?
            # Equation of the plane is r.n = a.n
            # Equation of the line is r = p + t(q-p) (where 0<t<1)
            # So the value of t for intersection is:
            #    t = (a.n - p.n) / ((q - p).n)
            dispDOTn = np.dot(nLocal, disp)
            t = (aDOTn[np.newaxis, :] - pDOTn) / dispDOTn[np.newaxis,:]
            intersects_plane = np.logical_and(t>0, t<1)
            intersectIds = intersects_plane.nonzero()

            # Now, put t back into the equation for the plane
            intersection = pGrid[intersectIds[0]] + np.outer(t[intersectIds][:,np.newaxis], disp)
            # Make relative to a
            intersection -= aLocal[intersectIds[1]]
            # Project to 2D by dropping the axis the triangle is most normal to.
            ignored_axis = (nLocal**2).argmax(axis=1)[intersectIds[1]]
            mask = np.ones_like(intersection, dtype=bool)
            mask[(np.arange(len(intersection)), ignored_axis)] = False
            AR = intersection[mask]
            AR.shape = (len(intersection), 2)
            AB = self.ab[triIds][intersectIds[1]][mask]
            AB.shape = AR.shape
            AC = self.ac[triIds][intersectIds[1]][mask]
            AC.shape = AR.shape
            # Is this within the triangle?
            # Find x & y s.t. x AB + y AC = AR
            # I.e.
            # /AB0 AC0\ /X\ = /AR0\
            # \AB1 AC1/ \Y/   \AR1/
            # If x > 0 and y > 0 and x+y < 1 it is inside.
            # Solutions are:
            # /X\ = ________1________ / AC1 -AC0\ /AR0\
            # \Y/   AB0.AC1 - AC0.AB1 \-AB1  AB0/ \AR1/
            det = (AB[:,0] * AC[:,1] - AC[:,0] * AB[:, 1]) 
            x = AC[:,1] * AR[:, 0] - AC[:,0] * AR[:,1]
            x /= det
            y = AB[:,0]*AR[:,1] - AB[:,1]*AR[:,0]
            y /= det
            inside_tri = np.logical_and(np.logical_and(x >= 0.0, y >= 0.0),
                                        x + y <= 1.0)
            # True for those hits where the direction p -> q is outward (for that triangle).
            hit_is_outward = dispDOTn[intersectIds[1][inside_tri]] > 0.0

            hit_pt_ids = intersectIds[0][inside_tri]
            hit_tri_loc_ids = intersectIds[1][inside_tri]
            hit_tri_ids = triIds[hit_tri_loc_ids]
            nHits = len(hit_pt_ids)
            total_hits += nHits
            
            for iHit in xrange(nHits):
                start_coord = np.array(np.unravel_index(hit_pt_ids[iHit], (size,size,size)))
                voxel_node = tri_node.GetNode(tri_node.levels, start_coord, create=True, relative=True)
                assert np.all(voxel_node.offset - tri_node.offset == start_coord)
                try:
                    intersections = voxel_node.intersections
                except AttributeError:
                    intersections = voxel_node.intersections = {}
                    pass

                # A hit is a tuple (t, outward_flag, tri_ID)
                new_hit = (t[hit_pt_ids[iHit], hit_tri_loc_ids[iHit]],
                             hit_is_outward[iHit],
                              hit_tri_ids[iHit])
                try:
                    hits = intersections[iDisp]
                    # We have previous hits for this link.
                    # Need to look for edge hits and merge them.
                    # We're just going to ignore the second if there is one.
                    merged = False
                    for iOldHit, h in enumerate(hits):
                        if abs(new_hit[0] - h[0]) < 1e-9:
                            # They are merge candidates
                            #old_id = h[2]
                            #new_id = hit_tri_ids[iHit]
                            if new_hit[1] == h[1]:
                                # Same sense, can merge these hits.
                                if h[0] > new_hit[0]:
                                    # New one closer, replace the old
                                    hits[iOldHit] = new_hit
                                merged = True
                                break
                            else:
                                #Let's hope this doesn't happen
                                pdb.set_trace()
                                pass
                            pass
                        continue
                    
                    if not merged:
                        hits.append(new_hit)
                except KeyError:
                    intersections[iDisp] = [new_hit]
                    pass
                continue
            continue
        return total_hits
    