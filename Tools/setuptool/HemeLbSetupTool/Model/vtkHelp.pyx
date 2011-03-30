cdef vtkOBBTree* FromPython(locator):
    cdef str empty, ptr, p, clsname
    empty, ptr, p, clsname = locator.__this__.split('_')
    assert clsname == 'vtkOBBTree'
    cdef unsigned int addr = int(ptr, 16)
    return <vtkOBBTree*>addr
