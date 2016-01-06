cimport cxdr

cdef class Unpacker:
    cdef cxdr.XDR x
    # def __cinit__(self, bytes buf)
    # def __dealloc__(self)
    cpdef double unpack_double(self)
    cpdef float unpack_float(self)
    cpdef int unpack_int(self)
    cpdef unsigned int unpack_uint(self)
