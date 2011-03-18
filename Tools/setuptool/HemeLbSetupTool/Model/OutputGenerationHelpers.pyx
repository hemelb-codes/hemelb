# cython: profile=True
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector

cimport numpy as np
import numpy as np

from HemeLbSetupTool.Model.Iolets import Inlet, Outlet
include "Flags.pxi"

# This ordering of the lattice vectors must be the same as in the HemeLB source.
cdef np.ndarray neighbours = np.array([[ 1, 0, 0],
                                       [-1, 0, 0],
                                       [ 0, 1, 0],
                                       [ 0,-1, 0],
                                       [ 0, 0, 1],
                                       [ 0, 0,-1],
                                       [ 1, 1, 1],
                                       [-1,-1,-1],
                                       [ 1, 1,-1],
                                       [-1,-1, 1],
                                       [ 1,-1, 1],
                                       [-1, 1,-1],
                                       [ 1,-1,-1],
                                       [-1, 1, 1]], dtype=np.int)
latticeVectorNorms = np.sqrt(np.sum(neighbours ** 2, axis= -1))
cdef np.ndarray _latticeVectorNorms = latticeVectorNorms
cdef np.ndarray laterNeighbourInds = np.array([0, 2, 4, 6, 8, 10, 12], dtype=np.uint)

# cdef struct Index:
#     unsigned int i, j, k
#     pass

ctypedef object Index
ctypedef object Domain

RealCode = np.float
ctypedef np.float_t Real
# ctypedef vector[np.float_t] RealVector

cdef class MacroBlock:
    cdef public Domain domain
    cdef public unsigned int size
    cdef public Index ijk
    cdef public np.ndarray sites
    
    def __init__(self, Domain domain, Index ijk, unsigned int size):
        self.domain = domain
        self.size = size
        self.ijk = ijk
        self.sites = np.empty((size, size, size), dtype=object)
        
        cdef Index globalSiteIjk
        
        for i from 0 <= i < size:
            for j from 0 <= j < size:
                for k from 0 <= k < size:
                    globalSiteIjk = (self.ijk[0] * size + i,
                                     self.ijk[1] * size + j,
                                     self.ijk[2] * size + k)
                    self.sites[i,j,k] = LatticeSite(self, globalSiteIjk)
        return

    cpdef IterSites(self):
        return self.sites.flat

    cpdef NdEnumerateSites(self):
        return np.ndenumerate(self.sites)
    
    cpdef GetLocalSite(self, int i, int j, int k):
        return self.sites[i, j, k]

    cpdef GetSite(self, int gi, int gj, int gk):
        cdef int i,j,k
        cdef int size = self.size
        i = gi - self.ijk[0]*size
        j = gj - self.ijk[1]*size
        k = gk - self.ijk[2]*size
        # Check if the coords belong to another block, i.e. any of
        # the local ones outside the range [0, self.size)
        if (i<0) or (j<0) or (k<0) or (i>=size) or (j>=size) or (k>=size):
            return self.domain.GetSite(gi, gj, gk)

        return self.GetLocalSite(i, j, k)
    
    cpdef CalcPositionFromIndex(self, Index index):
        return self.domain.CalcPositionFromIndex(index)

cdef class LatticeSite:
    cdef public MacroBlock block
    cdef public Index ijk
    cdef vector[Real]* _Position
    
    cdef public bint IsFluid
    cdef public object IsEdge
    cdef public object Iolet
    cdef public object BoundaryId
        
    # Attributes that will be updated by Profile.ClassifySite
    cdef public np.ndarray BoundaryNormal
    cdef public Real BoundaryDistance
    
    cdef public np.ndarray WallNormal
    cdef public Real WallDistance
    
    cdef public np.ndarray CutDistances
    cdef public np.ndarray CutCellIds
    
    def __cinit__(self, MacroBlock block, Index ijk):
        self._Position = new vector[Real](3)
        self.block = block
        self.ijk = ijk
        pos = block.CalcPositionFromIndex(ijk)
        cdef Real val
        # cdef vector[Real] pv = deref(self._Position)
        for i from 0 <= i < 3:
            # val = pos[i]
            deref(self._Position)[i] = pos[i]
            

        self.IsFluid = False
        self.IsEdge = None
        self.Iolet = None
        self.BoundaryId = None

        # Attributes that will be updated by Profile.ClassifySite
#        self.Type = None
        self.BoundaryNormal = np.zeros(3)
        self.BoundaryDistance = float('inf')

        self.WallNormal = np.zeros(3)
        self.WallDistance = float('inf')

        self.CutDistances = 1. / np.zeros(len(neighbours), dtype=RealCode)
        self.CutCellIds = -np.ones(len(neighbours), dtype=int)
        return
    
    property Position:
        def __get__(self):
            cdef np.ndarray ans = np.zeros(3, dtype=RealCode)
            for i from 0 <= i < 3:
                ans[i] = deref(self._Position)[i]
            return ans
        def __set__(self, np.ndarray val):
            for i from 0 <= i < 3:
                deref(self._Position)[i] = val[i]
    
    property Type:
        def __get__(self):
            if self.IsFluid:
                if isinstance(self.Iolet, Inlet):
                    return INLET_TYPE
                elif isinstance(self.Iolet, Outlet):
                    return OUTLET_TYPE
                else:
                    return FLUID_TYPE
            else:
                return SOLID_TYPE
    
    property Config:
        def __get__(self):
            cfg = self.Type
            if not self.IsFluid:
                return self.Type

            # Fluid sites now
            if self.Type == FLUID_TYPE and not self.IsEdge:
                # Simple fluid sites
                return self.Type

            # A complex one

            if self.IsEdge:
                # Bit fiddle the boundary config. See comment below
                # by BOUNDARY_CONFIG_MASK for definition.
                boundary = 0
                for i, neigh in self.EnumerateNeighbours():
                    if not neigh.IsFluid:
                        # If the lattice vector is cut, set the flag
                        boundary |= 1 << i
                        pass
                    continue
                # Shift the boundary bit field to the appropriate
                # place and set these bits in the type
                cfg |= boundary << BOUNDARY_CONFIG_SHIFT

                if boundary:
                    # Set this bit if we've hit any solid sites
                    cfg |= PRESSURE_EDGE_MASK
                    pass
                pass

            if self.Type != FLUID_TYPE:
                # It must be an inlet or outlet
                # Shift the index left and set the bits
                cfg |= self.BoundaryId << BOUNDARY_ID_SHIFT
                pass

            return cfg

    def EnumerateNeighbours(self):
        """Create an iterator that enumerates our neighbours.
        """
        return SiteNeighbourEnumerator(self)

    def EnumerateLaterNeighbours(self):
        """Create an iterator that enumerates our neighbours who're later in the array.
        """
        return SiteLaterNeighbourEnumerator(self)
    pass

cdef class SiteNeighbourEnumeratorBase:
    cdef object domain
    cdef object shape
    cdef np.int_t i
    cdef np.int_t max_i
    cdef LatticeSite site
    
    def __cinit__(self, LatticeSite site):
        self.site = site
        self.domain = site.block.domain
        self.shape = self.domain.BlockCounts * self.domain.BlockSize
        self.i = -1
        return
    
    cdef np.ndarray GetVector(SiteNeighbourEnumeratorBase self):
        return neighbours[self.i]
    
    def __iter__(self): return self
    
    def __next__(self):
        cdef np.ndarray latticeVec
        cdef SiteNeighbourEnumeratorBase slf = self
        while True:
            slf.i += 1
            if slf.i >= slf.max_i:
                raise StopIteration
            
            latticeVec = slf.GetVector()
            ind = slf.site.ijk + latticeVec
            # Skip any out of range
            if np.any(ind < 0) or np.any(ind >= self.shape):
                continue
            break
        
        return self.i, self.site.block.GetSite(ind[0], ind[1], ind[2])
    pass

cdef class SiteNeighbourEnumerator(SiteNeighbourEnumeratorBase):
    def __cinit__(self, LatticeSite site):
        self.max_i = 14
        return
    pass

cdef class SiteLaterNeighbourEnumerator(SiteNeighbourEnumeratorBase):
    def __cinit__(self, LatticeSite site):
        self.max_i = 7
        return
    cdef np.ndarray GetVector(self):
        return neighbours[laterNeighbourInds[self.i]]
    pass

