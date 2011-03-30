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

cdef vtkOBBTree* FromPython(locator)
