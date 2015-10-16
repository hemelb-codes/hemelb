import numpy as np
from itertools import chain        

from CGAL.CGAL_Kernel import Point_3, Segment_3
from CGAL.CGAL_AABB_tree import AABB_tree_Polyhedron_3_Facet_handle
from CGAL.CGAL_Polyhedron_3 import Polyhedron_3

from Neighbours import neighbours, inverses, norm2

class Face(object):
    def __init__(self, size, norm):
        self.iter_dims = []
        for i, n in enumerate(norm):
            if n == 0:
                self.iter_dims.append(i)
            else:
                self.norm_dim = i
                
        self.normal = np.array(norm, dtype=int)
        self.size = size
        self.name = {+1: '+', -1: '-'}[norm[self.norm_dim]] 
        for iter_dim in self.iter_dims:
            self.name += {0: 'x', 1: 'y', 2: 'z'}[iter_dim]
        
        self.coord = {+1: -1, -1: size}[norm[self.norm_dim]]
        return
     
    def dot(self, vec):
        return np.sum(self.normal * vec)
            
    def __iter__(self):
        return self.Iter()
    
    def Iter(self, *prev_faces):
        limits =[ [-1, self.size+1],
                  [-1, self.size+1],
                  [-1, self.size+1] ]
        for f in prev_faces:
            i, val = {
                      -1: (0, 0),
                      self.size: (1, self.size)
                      }[f.coord]
            limits[f.norm_dim][i] = val
            continue
        limits[self.norm_dim] = (self.coord, self.coord+1)
        
        for i in xrange(*limits[0]):
            for j in xrange(*limits[1]):
                for k in xrange(*limits[2]):
                    yield np.array((i,j,k), dtype=int)

class FaceGroup(object):
    def __init__(self, *faces):
        self.faces = faces
    def __iter__(self):
        face_iterators = [self.faces[i].Iter(*self.faces[:i])
                          for i in xrange(len(self.faces))]
        return chain(*face_iterators)
    pass

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
        
        pts = self.CGALpoints = np.empty(len(points), dtype=object) 
        for i, p in enumerate(points):
            pts[i] = Point_3(*p)
        # Supply number of vertices, half-edges and facets
        surf = self.CGAL_surface = Polyhedron_3(len(points),
                                                3*len(triangles),
                                                len(triangles))
        
        for i, tri in enumerate(triangles):
            half_edge = surf.make_triangle(pts[tri[0]],
                                           pts[tri[1]],
                                           pts[tri[2]])
            facet = half_edge.facet()
            facet.set_id(i)
            
        self.aabb_tree = AABB_tree_Polyhedron_3_Facet_handle(self.CGAL_surface.facets())
        return
    
    def __call__(self, tri_node):
        size = 2**tri_node.levels
        faces = []
        for d in xrange(3):
            for pm_one in (+1, -1):
                n = np.zeros(3, dtype=int)
                n[d] = pm_one
                faces.append(Face(size, n))
                
        for i_vec, vec in enumerate(neighbours):
            if inverses[i_vec] < i_vec:
                continue
            fg = FaceGroup(*[f for f in faces if f.dot(vec) > 0])
            for ijk in fg:
                i = 1
                node_coord = ijk + i*vec
                while np.all(np.logical_and(node_coord>=0, node_coord<size)): 
                    i += 1
                    node_coord = ijk + i*vec
                    continue
                if i > 1:
                    start = ijk + tri_node.offset
                    end = node_coord + tri_node.offset
                    self._DoRay(tri_node, i_vec, start, end)
                    pass
                continue
            continue
        return
    
    def _DoRay(self, tri_node, i_vec, startIjk, endIjk):
        start = Point_3(*startIjk) #*[float(x) for x in (startIjk + tri_node.offset)])
        end = Point_3(*endIjk)
        ray = Segment_3(start, end)
        intersections = []
        self.aabb_tree.all_intersections(ray, intersections)
        vec = neighbours[i_vec]
        for pair in intersections:
            CGAL_ip = pair.first.get_Point_3()
            r_hit = np.array((CGAL_ip.x(), CGAL_ip.y(), CGAL_ip.z()))
            dr = r_hit - startIjk
            alongness = np.sum(dr*vec) / norm2[i_vec]
            
            vox1 = int(alongness)*vec + startIjk
            t1 = alongness - int(alongness)
            triId = pair.second.id()
            out12 = np.sum(self.normals[triId] * vec) >= 0.0
            # A hit is a tuple (t, outward_flag, tri_ID)
            # Need to deal with the case where vox is outside this node's ROI
            if not np.all(vox1 == startIjk):
                self.add_hit(tri_node, vox1, i_vec,
                        (t1, out12, triId))
            vox2 = vox1 + vec
            # And the same for vox2
            if not np.all(vox2 == endIjk):
                self.add_hit(tri_node, vox2, inverses[i_vec],
                        (1.0 - t1, not out12, triId))
            
            continue
        return
    
    @staticmethod
    def add_hit(tri_node, ijk, i_vec, new_hit):
        vox = tri_node.GetNode(0, ijk, create=True)
        try:
            intersections = vox.intersections
        except AttributeError:
            intersections = vox.intersections = {}
            pass
        # A hit is a tuple (t, outward_flag, tri_ID)
        try:
            hits = intersections[i_vec]
            hits.append(new_hit)
        except KeyError:
            intersections[i_vec] = [new_hit]
            pass
        return
    pass
    