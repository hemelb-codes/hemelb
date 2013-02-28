# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 
cimport numpy as np
import numpy as np
import weakref

cdef public enum:
    DIRECTIONS = 26
    
cdef public enum:
    SOLID = 0
    FLUID = 1
    
cdef public enum:
    NO_INTERSECTION = 0
    WALL_INTERSECTION = 1
    INLET_INTERSECTION = 2
    OUTLET_INTERSECTION = 3

cdef public enum:
    NO_IOLET = -1

cdef public enum:
    WALL_NORMAL_NOT_AVAILABLE = 0
    WALL_NORMAL_AVAILABLE = 1

cdef class BaseSite:
    def __init__(self, block, np.ndarray[np.int_t] sgIdx):
        self.GetBlock = weakref.ref(block)
        self.Index = sgIdx
        self.IntersectionType = None
        self.IntersectionDistance = None
        self.IOletIndex = None
        self.Type = SOLID
        self.WallNormalAvailable = False
        self.WallNormal = None

    cdef bint _IsFluid(self):
        return self.Type == FLUID
    property IsFluid:
        def __get__(self):
            return self._IsFluid()
        
    cdef bint _IsSolid(self):
        return self.Type == SOLID
    property IsSolid:
        def __get__(self):
            return self._IsSolid()
    
    cdef bint _IsEdge(self):
        return self.IsFluid and np.any(self.IntersectionType == WALL_INTERSECTION)
    property IsEdge:
        def __get__(self):
            return self._IsEdge()
    
    cdef np.ndarray _Position(self):
        dom = self.GetBlock().GetDomain()
        return dom.Origin + dom.VoxelSize * self.Index
    property Position:
        def __get__(self):
            return self._Position()

    def __getstate__(self):
        picdic = {
            'Type': self.Type,
            'Index': self.Index,
            'IntersectionType': self.IntersectionType,
            'IntersectionDistance': self.IntersectionDistance,
            'IOletIndex': self.IOletIndex,
            'Block': self.GetBlock(),
            'WallNormalAvailable': self.WallNormalAvailable,
            'WallNormal': self.WallNormal
            }
        return picdic

    def __setstate__(self, picdic):
        self.Type = picdic['Type']
        self.Index = picdic['Index']
        self.IntersectionType = picdic['IntersectionType']
        self.IntersectionDistance = picdic['IntersectionDistance']
        self.IOletIndex = picdic['IOletIndex']
        self.GetBlock = weakref.ref(picdic['Block'])
        self.WallNormalAvailable =  picdic['WallNormalAvailable']
        self.WallNormal = picdic['WallNormal']

    cpdef LoadFrom(self, xdr.Unpacker loader):
        cdef np.ndarray[np.uint_t] itype
        cdef np.ndarray[np.float32_t] idist
        cdef np.ndarray[np.int_t] ioind
        cdef int i
        
        self.Type = loader.unpack_uint()
        # Solid and simple fluid, we are done loading
        if self.IsFluid == FLUID:
            # Initialise arrays
            itype = self.IntersectionType = np.empty(DIRECTIONS, dtype=np.uint)
            itype[:] = NO_INTERSECTION
            idist = self.IntersectionDistance = np.zeros(DIRECTIONS, dtype=np.float32)
            ioind = self.IOletIndex = np.empty(DIRECTIONS, dtype=np.int)
            ioind[:] = NO_IOLET
            
            for i in xrange(DIRECTIONS):
                itype[i] = loader.unpack_uint()
                if itype[i] == INLET_INTERSECTION or itype[i] == OUTLET_INTERSECTION:
                    ioind[i] = loader.unpack_uint()
                if itype[i] != NO_INTERSECTION:
                    idist[i] = loader.unpack_float()
            
            self.WallNormalAvailable = (loader.unpack_uint() == WALL_NORMAL_AVAILABLE)
            if self.WallNormalAvailable:
                self.WallNormal = np.array([loader.unpack_float(), loader.unpack_float(), loader.unpack_float()], dtype=np.float32)
