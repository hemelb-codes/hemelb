import numpy as np
    
from vtkHelp cimport vtkPoints_New, vtkIdList_New

cdef class PointCellIdPair:
    def __cinit__(self):
        self.p = None

cdef class HitList(object):
    """Object that finds the intersections given by a vtkLocator the
    line joining a start and end LatticeSite. It then gives a simple
    iterator-like method (next) for getting these points, which when
    it finishes, returns None rather than raise StopIteration.
    """
    
    def __cinit__(self):
        self.HitPoints = vtkPoints_New()
        self.HitCellIds = vtkIdList_New()
        
        self.Point = np.zeros(3, dtype=np.float)

    cdef Init(self, vtkOBBTree* locator, np.ndarray start, np.ndarray end):
        """locator - vtkLocator to use to test for intersections.
        start - starting LatticeSite
        end - ending LatticeSite
        """
        # Find intersections
        locator.IntersectWithLine(<double*>start.data,
                              <double*>end.data,
                              self.HitPoints, self.HitCellIds)
        # Cache the number of them
        self.len = self.HitPoints.GetNumberOfPoints()
        # Track where we are in the iteration
        self.i = 0
        
    
    def __dealloc__(self):
        self.HitPoints.Delete()
        self.HitCellIds.Delete()

    cdef PointCellIdPair cNext(self):
        cdef PointCellIdPair ans = PointCellIdPair()
        cdef int i = self.i
        if i < self.len:
            self.i += 1
            
            self.HitPoints.GetPoint(i, <double*>self.Point.data)
            ans.p = self.Point
            ans.id = self.HitCellIds.GetId(i)
            return ans
        else:
            raise StopIteration

    def __next__(self):
        """Move our internal internal counter on one and return the
        new value, or None in the case that it doesn't exist.
        """
        ans = self.cNext()
        return ans.p, ans.id
        
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
