# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#
cimport cxdr

cdef class Unpacker:
    cdef cxdr.XDR x
    # def __cinit__(self, bytes buf)
    # def __dealloc__(self)
    cpdef double unpack_double(self)
    cpdef float unpack_float(self)
    cpdef int unpack_int(self)
    cpdef unsigned int unpack_uint(self)
