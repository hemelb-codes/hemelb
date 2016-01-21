
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
cimport cxdr

cdef class Unpacker:
    cdef cxdr.XDR x
    # def __cinit__(self, bytes buf)
    # def __dealloc__(self)
    cpdef double unpack_double(self)
    cpdef float unpack_float(self)
    cpdef int unpack_int(self)
    cpdef unsigned int unpack_uint(self)
