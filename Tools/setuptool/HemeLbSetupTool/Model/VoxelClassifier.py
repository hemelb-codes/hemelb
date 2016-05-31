#
# Copyright (C) University College London, 2007-2012, all rights reserved.
#
# This file is part of HemeLB and is CONFIDENTIAL. You may not work
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#

import numpy as np
from Oct import NodeError
import HemeOct
from . import Lattice
from pytest import set_trace

class Cut(object):
    def __init__(self, direction, cutdist, triId):
        self.links = {}
        self.dist = {direction: cutdist,
                     Lattice.Opposites[direction]: 1-cutdist}
        self.triId = triId
        return
    pass

class PlaceholderVoxel(object):
    def __init__(self, offset):
        self.offset = offset
        self.links = {}
    
class LocalPlaceholderVoxel(PlaceholderVoxel):
    pass

class NonLocalPlaceholderVoxel(PlaceholderVoxel):
    pass

class DLL(object):
    """Doubly linked list along links.
    """
    def __init__(self, direction, start):
        self.head = start
        self.tail = start
        self.i = direction
        self.iOpp = Lattice.Opposites[direction]
        
        if not hasattr(self.head, "links"):
            self.head.links = {}
            
        assert not self.head.links.has_key(self.i)    
        self.head.links[self.i] = None
        
    def append(self, obj):
        if not hasattr(obj, "links"):
            obj.links = {}
            
        obj.links[self.i] = None
        obj.links[self.iOpp] = self.tail
        
        self.tail.links[self.i] = obj
        self.tail = obj
    
    def __len__(self):
        n = 0
        ptr = self.head
        while ptr is not None:
            n += 1
            ptr = ptr.links[self.i]
        return n
    
class VoxelClassifier(object):
        
    def __init__(self, points, triangles, normals, labels):
        self.Points = points
        self.Triangles = triangles
        self.Normals = normals
        self.Labels = labels
        return
    
    def __call__(self, tree, tri_level):
        ans = HemeOct.Tree(tree.levels)
        for node in tree.IterDepthFirst(tri_level, tri_level):
            region = self.DoRegion(node)
            if region.children is not None or np.any(region.children) is not None:
                ans.root.SetNode(tri_level, node.offset, region)
                pass
            continue
        return ans
            
    def IntersectLinkWithTriangle(self, coord, displacement, iTri):
        tri_pt_ids = self.Triangles[iTri]
        tri_points = self.Points[tri_pt_ids]
        # Work relative to the voxel coord
        tri_points -= coord
        norm = self.Normals[iTri]
        # assume that the line intersects the triangle's plane at a fraction t 
        # along the displacement vector
        t = np.dot(tri_points[0], norm) / np.dot(displacement, norm)
        # To intersect, t must be in (0,1)
        if t < 0 or t > 1:
            return np.inf
        # so the intersection is on our line segment, is it in the triangle?
        # first find the point
        r = t * displacement
        # get coords relative to one point of the tri
        v10 = tri_points[1] - tri_points[0]
        v20 = tri_points[2] - tri_points[0]
        vr0 =  r - tri_points[0]

        v10_v10 = np.sum(v10**2)
        v20_v20 = np.sum(v20**2)
        v10_v20 = np.dot(v10, v20)
        
        denom = v10_v10 * v20_v20 - v10_v20**2
        
        vr0_v10 = np.dot(vr0, v10)
        vr0_v20 = np.dot(vr0, v20)
        # compute the barycentric coords of the intersection
        u = (v20_v20 * vr0_v10 - v10_v20 * vr0_v20) / denom
        v = (v10_v10 * vr0_v20 - v10_v20 * vr0_v10) / denom
        
        if u>=0 and v>=0 and u+v<=1:
            return t
        
        return np.inf
    
    def DoRegion(self, region_node):
        """This accepts a subtree of voxels that represent the surface and have 
        data about which triangles may intersect their links.
        
        It examines all links from this set looking for intersections.
        """
        region_extras = {}
        other_extras = {}
        new_region = HemeOct.Node(region_node.levels, region_node.offset)
        for voxel in region_node.IterDepthFirst(0,0):
            if not hasattr(voxel, 'links'):
                links = voxel.links = {}
                pass
            
            for i, disp in enumerate(Lattice.Neighbours):
                if voxel.links.has_key(i):
                    # We've done this link from the other direction, skip
                    continue
                neighIdx = voxel.offset + disp
                
                if np.all(voxel.offset == (9,13,22)) and np.all(neighIdx == (8,14,21)):
                    set_trace()
                if np.all(neighIdx == (9,13,22)) and np.all(voxel.offset == (8,14,21)):
                    set_trace()
                
                # compute all cuts
                cuts = []
                for tri in voxel.triIds:
                    cut_dist = self.IntersectLinkWithTriangle(voxel.offset, disp, tri)
                    if cut_dist <= 1.0:
                        cuts.append(Cut(i, cut_dist, tri))
                
                # Sort this list
                cuts.sort(key=lambda c: c.dist[i])
                cut_dll = DLL(i, voxel)
                # Add them to the vox-vox linked list
                for cut in cuts:
                    cut_dll.append(cut)
                    continue
                
                try:
                    neigh = region_node.GetNode(0, neighIdx)
                except IndexError:
                    # node out of range
                    neigh = NonLocalPlaceholderVoxel(neighIdx)
                    other_extras[tuple(neighIdx)] = neigh
                    pass
                except NodeError:
                    # Node doesn't exist, but could
                    neigh = LocalPlaceholderVoxel(neighIdx)
                    region_extras[tuple(neighIdx)] = neigh
                    pass
 
                cut_dll.append(neigh)
                continue # End loop over links
            
            # Look at all my cuts, counting the number of innies and outies
            nIn = 0
            nOut = 0
            for i, cut in voxel.links.iteritems():
                if isinstance(cut, Cut):
                    norm = self.Normals[cut.triId]
                    if np.dot(Lattice.Neighbours[i], norm) > 0:
                        nOut += 1
                    else:
                        nIn += 1
            if nOut > 0:
                if nIn > 0:
                    # Oh dear - we have an inconsistency
                    raise ValueError("Displacement dot surface normal gives inconsistent results for point at " + str(voxel.offset))
                # We are inside
                new_voxel = new_region.GetNode(0, voxel.offset, create=True)
                new_voxel.fluid_count = 1
                
                new_voxel.cut_tri = np.zeros(26, dtype=int)
                new_voxel.cut_tri[:] = -1
                
                new_voxel.cut_dist = np.zeros(26, dtype=float)
                for i, cut in voxel.links.iteritems():
                    if isinstance(cut, Cut):
                        new_voxel.cut_tri[i] = cut.triId
                        new_voxel.cut_dist[i] = cut.dist[i]
                        
            elif nIn > 0:
                # We are outside, don't add to output
                pass
            else:
                # No hits at all don't add this to the output
                pass
            
            continue # End of iter over voxels that existed at start
        
        # Now deal with the extra voxels
        return new_region
    pass
