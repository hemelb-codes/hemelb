# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 
cimport numpy as np
from hemeTools.utils cimport xdr

cdef class BaseSite:
    cdef public:
        int Type
        np.ndarray Index
        np.ndarray IntersectionType
        np.ndarray IntersectionDistance
        np.ndarray IoletIndex
        object GetBlock

    cdef bint _IsFluid(self)
    cdef bint _IsSolid(self)
    cdef bint _IsEdge(self)
    cdef np.ndarray _Position(self)
    cpdef LoadFrom(self, xdr.Unpacker loader)
    
    pass
