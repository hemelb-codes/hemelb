cimport numpy as np
from vtkHelp cimport vtkOBBTree, vtkPoints, vtkIdList

cdef class PointCellIdPair:
    cdef:
        np.ndarray p
        int id

cdef class HitList(object):
    cdef:
        vtkPoints* HitPoints
        vtkIdList* HitCellIds
        np.ndarray Point
        int len
        int i
        
    cdef Init(self, vtkOBBTree* locator, np.ndarray start, np.ndarray end)
    cdef PointCellIdPair cNext(self)
