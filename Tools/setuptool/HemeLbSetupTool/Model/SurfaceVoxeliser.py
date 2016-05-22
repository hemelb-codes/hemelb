import numpy as np
import itertools
import multiprocessing
import HemeOct
from pytest import set_trace

class SurfaceVoxeliser(object):
    """This class creates an octree with leaf nodes that make a minimal, 
    26-separable voxelisation of the surface represented by the triangle mesh, 
    following Huang et al (http://dx.doi.org/10.1109/SVV.1998.729593)
    
    The leaf nodes have added to them a set listing the triangles that they 
    could intersect (triIds)
    """
    
    def __init__(self, points, triangles, normals, tree, tri_level):
        """
        points = a numpy array of doubles with shape == (nPoints, 3) specifying
        the vertices.
         
        triangles = a numpy array of ints with shape == (nTris, 3) specifiying
        the IDs of the verticies comprising the triangle.
        
        normals = a numpy array of doubles with shape == (nPoints, 3) specifying 
        the normals
        
        tree = the octree with the triangles sorted onto nodes at a given level
        
        tri_level = the level of the tree where the triIds attribute lives
        
        """
        self.Points = points
        self.Triangles = triangles
        self.Normals = normals
                
        self.Tree = tree
        self.TriLevel = tri_level
        
        # This is in units of the voxel size
        self.Rc = np.sqrt(3.0) / 2 
        
    def AABB(self, iTri):
        """Return the axis-aligned bounding box for a triangle."""
        ids = self.Triangles[iTri]
        points = self.Points[ids]
        lo = points.min(axis=0)
        hi = points.max(axis=0)
        return lo, hi
    
    def FilterPoint(self, iPt, voxels, inside_mask):
        """Create all voxels lying within Rc = sqrt(3)/2 of the point"""
        
        p = self.Points[iPt]
        dr = voxels - p
        dr2 = np.sum(dr**2,axis=-1)
        
        inside_mask |= dr2 < 0.75
        return
    
    def FilterEdge(self, iPt, jPt, voxels, inside_mask):
        """
        x .              . b
        
        
               .  p       
        
        a .
        
        Our point of interest is x
        line is a -> b
        it's equation is a + lambda n
        where 0 <= lambda <= 1 
        and n = b - a
        p is the closest point on the line to x
        define r = x - a
        
        Since xp is perpendicular to ab:
            p.n = x.n
        hence
            lambda = (x - a).n / n.n
        and
            p = a + n*(r.n) / (n.n)
             
        if rho = |x - p|, Pythagoras tells us:
            r.r = rho^2 + (r.n)^2 / (n.n)
        """
        a = self.Points[iPt]
        b = self.Points[jPt]
        
        n = b - a
        n2 = np.dot(n, n)
        
        r = voxels - a
        r_n = np.dot(r, n)
        lambdas =  r_n / n2
        
        r2 = np.sum(r**2, axis=-1)
        rho2 = r2 - r_n**2 / n2 
        
        # Inside points have 0 <= lamdba <= 1 and rho^2 < Rc^2
        inside_mask |= (0 <= lambdas) & (lambdas <= 1) & (rho2 <= 0.75)
        return
    
    def FilterPlane(self, iTri, voxels, inside_mask):
        """Mark as inside all points within the triangular prism defined by the
        following 2 planes (see Huang Fig. 12)
        """
        ids = self.Triangles[iTri]
        norm = self.Normals[iTri]
        points = self.Points[ids]
        
        # t_26 is defined in Huang figure 10
        t26 = np.sum(np.abs(norm) /2)

        upper_points = points + t26*norm
        lower_points = points - t26*norm

        # We want the *inward* normals        
        plane_normals = np.zeros((2,3), dtype=float)
        plane_points = np.zeros((2,3), dtype=float)
        
        # Top and bottom are easy
        plane_normals[0] = norm
        plane_points[0] = lower_points[0]
        
        plane_normals[1] = -norm
        plane_points[1] = upper_points[0]
       
        plane_offsets = np.sum(plane_normals * plane_points, axis=-1)
        
        # greater than because the normals are inwards
        inside_mask |= np.all(np.dot(voxels, plane_normals.transpose()) > plane_offsets,
                              axis=1)
        return

    def FilterTriangle(self, iTri, voxels, inside_mask):
        """Mark as inside all points within the triangular prism defined by the
        following 5 planes (see Huang Fig. 12)
        """
        ids = self.Triangles[iTri]
        norm = self.Normals[iTri]
        points = self.Points[ids]
        
        # t_26 is defined in Huang figure 10
        t26 = np.sum(np.abs(norm) /2)

        upper_points = points + t26*norm
        lower_points = points - t26*norm

        # We want the *inward* normals        
        plane_normals = np.zeros((5,3), dtype=float)
        plane_points = np.zeros((5,3), dtype=float)
        
        # Top and bottom are easy
        plane_normals[0] = norm
        plane_points[0] = lower_points[0]
        
        plane_normals[1] = -norm
        plane_points[1] = upper_points[0]
        
        # Think about the edges now
        # Will need cross products to get their normals. Our normal can be used
        # as one of the inputs.
        v10 = points[1] - points[0]
        v21 = points[2] - points[1]
        v02 = points[0] - points[2]
        
        plane_normals[2] = np.cross(norm, v10)
        plane_points[2] = upper_points[0]
        
        plane_normals[3] = np.cross(norm, v21)
        plane_points[3] = upper_points[1]
        
        plane_normals[4] = np.cross(norm, v02)
        plane_points[4] = upper_points[2]
        
        plane_offsets = np.sum(plane_normals * plane_points, axis=-1)
        
        # greater than because the normals are inwards
        inside_mask |= np.all(np.dot(voxels, plane_normals.transpose()) > plane_offsets,
                              axis=1)
        return
    
    def DoTriangle(self, iTri):
        """Add the voxel nodes for a single triangle to the tree
        """
        norm = self.Normals[iTri]
        ids = self.Triangles[iTri]
        points = self.Points[ids]
        
        lo, hi = self.AABB(iTri)
        # voxels that could in principle intersect
        vlo = lo.astype(int)
        vhi = hi.astype(int) + 2
        vshape = vhi - vlo
        vsize = np.prod(vshape)
        # +2 for use as upper bound in range statement
        
        voxels = np.mgrid[vlo[0]:vhi[0],
                          vlo[1]:vhi[1],
                          vlo[2]:vhi[2]].reshape((3, vsize)).transpose()
        
        inside_mask = np.zeros(vsize, dtype=bool)
        
        # Now apply the three tests:
        # 1) within point sphere?
        # 2) within edge cylinder?
        # 3) within the planes?
        
        for i in xrange(3):
            iPt = ids[i]
            jPt = ids[(i + 1) % 3]
            self.FilterPoint(iPt, voxels, inside_mask)
            self.FilterEdge(iPt, jPt, voxels, inside_mask)
                    
        self.FilterTriangle(iTri, voxels, inside_mask)
        
        for vox in voxels[inside_mask.nonzero()]:
            node = self.Tree.GetNode(0, vox, create=True)
            try:
                node.triIds.add(iTri)
            except AttributeError:
                node.triIds = set((iTri,))
                
        
        return
    
    def Execute(self):
        for node in self.Tree.IterDepthFirst(self.TriLevel, self.TriLevel):
            ids = node.triIds
            for iTri in ids:
                self.DoTriangle(iTri)
                