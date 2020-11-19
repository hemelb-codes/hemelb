# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
cimport numpy as np
from hemeTools.utils cimport xdr

cdef class BaseSite:
    cdef public:
        # BEWARE: Adding attributes to this class requires programming
        # their pickling/unpickling in __getstate__ and __setstate__
        int Type
        np.ndarray Index
        np.ndarray IntersectionType
        np.ndarray IntersectionDistance
        np.ndarray IOletIndex
        bint WallNormalAvailable
        np.ndarray WallNormal
        object GetBlock

    cdef bint _IsFluid(self)
    cdef bint _IsSolid(self)
    cdef bint _IsEdge(self)
    cdef np.ndarray _Position(self)
    cpdef LoadFrom(self, xdr.Unpacker loader)
    
    pass
