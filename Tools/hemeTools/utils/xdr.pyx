# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#
cimport cxdr

cdef class Unpacker:

    def __cinit__(self, bytes buf):
        # cdef bytes cBuf = buf.encode()
        cxdr.xdrmem_create(&self.x, buf, len(buf), cxdr.XDR_DECODE)

    def __dealloc__(self):
        cxdr.xdr_destroy(&self.x)

    cpdef double unpack_double(self):
        cdef double ans
        cxdr.xdr_double(&self.x, &ans)
        return ans

    cpdef float unpack_float(self):
        cdef float ans
        cxdr.xdr_float(&self.x, &ans)
        return ans
    
    cpdef int unpack_int(self):
        cdef int ans
        cxdr.xdr_int(&self.x, &ans)
        return ans
    
    cpdef unsigned int unpack_uint(self):
        cdef unsigned int ans
        cxdr.xdr_u_int(&self.x, &ans)
        return ans
    
    
