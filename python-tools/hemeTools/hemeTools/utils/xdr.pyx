# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cimport cxdr
from libc.stdint cimport int32_t, uint32_t

cdef class Unpacker:
    def __cinit__(self, bytes buf):
        self._buf = buf
        self._pos = 0

    cpdef double unpack_double(self):
        cdef double ans
        cxdr.xdr_deserialise_double(&ans, &self._buf[self._pos])
        self._pos += 8
        return ans

    cpdef float unpack_float(self):
        cdef float ans
        cxdr.xdr_deserialise_float(&ans, &self._buf[self._pos])
        self._pos += 4
        return ans
    
    cpdef int unpack_int(self):
        cdef int32_t ans
        cxdr.xdr_deserialise_int32(&ans, &self._buf[self._pos])
        self._pos += 4
        return ans
    
    cpdef unsigned int unpack_uint(self):
        cdef uint32_t ans
        cxdr.xdr_deserialise_uint32(&ans, &self._buf[self._pos])
        self._pos += 4
        return ans
