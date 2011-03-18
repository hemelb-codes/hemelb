import numpy as np
cimport numpy as np

cdef extern from "vtkIdList.h":
    cdef cppclass vtkIdList:
        int GetId (int i)
        void Delete()
    cdef vtkIdList* vtkIdList_New "vtkIdList::New"()
    
cdef extern from "vtkPoints.h":
    cdef cppclass vtkPoints:
        int GetNumberOfPoints()
        void GetPoint (int id, double* x)
        void Delete()
    cdef vtkPoints* vtkPoints_New "vtkPoints::New"()
    
cdef extern from "vtkOBBTree.h":
    cdef cppclass vtkOBBTree:
        int IntersectWithLine(double* a0, double* a1, vtkPoints *points, vtkIdList *cellIds)
        
    cdef vtkOBBTree* vtkOBBTree_New "vtkOBBTree::New"()
    
    
cdef class HitList(object):
    """Object that finds the intersections given by a vtkLocator the
    line joining a start and end LatticeSite. It then gives a simple
    iterator-like method (next) for getting these points, which when
    it finishes, returns None rather than raise StopIteration.
    """
    cdef:
        vtkPoints* HitPoints
        vtkIdList* HitCellIds
        np.ndarray Point
        int len
        int i
        
    def __cinit__(self, locator, np.ndarray start, np.ndarray end):
        """locator - vtkLocator to use to test for intersections.
        start - starting LatticeSite
        end - ending LatticeSite
        """
        self.HitPoints = vtkPoints_New()
        self.HitCellIds = vtkIdList_New()
        
        self.Point = np.zeros(3, dtype=np.float)
        
        empty, ptr, p, clsname = locator.__this__.split('_')
        assert clsname == 'vtkOBBTree'
        cdef unsigned int addr = int(ptr, 16)
        cdef vtkOBBTree* loc = <vtkOBBTree*>addr
        
        # Find intersections
        loc.IntersectWithLine(<double*>start.data,
                              <double*>end.data,
                              self.HitPoints, self.HitCellIds)
        # Cache the number of them
        self.len = self.HitPoints.GetNumberOfPoints()
        # Track where we are in the iteration
        self.i = 0
        
    
    def __dealloc__(self):
        self.HitPoints.Delete()
        self.HitCellIds.Delete()
        
    def __next__(self):
        """Move our internal internal counter on one and return the
        new value, or None in the case that it doesn't exist.
        """
        i = self.i
        if i < self.len:
            self.i += 1
            return self[i]
        else:
            raise StopIteration

    def __iter__(self): return self
    
    def __len__(self):
        return self.len

    def __getitem__(self, int i):
        if i < 0: raise IndexError
        if i >= self.len: raise IndexError

        self.HitPoints.GetPoint(i, <double*>self.Point.data)

        return self.Point, self.HitCellIds.GetId(i)
    pass

# def IterHitsForLink(locator, start, end):
#         """Given a pair of sites, yield the itersections one at a
#         time.
#         """
#         # Intersect against the STL surface
#         stlHits = StlHitList(self.Locator, start, end)

#         # Now start yielding intersections
#         topStl = stlHits.next()

#         while topStl is not None:
#             yield topStl
#             topStl = stlHits.next()
#             continue

#         return
