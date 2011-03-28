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

RealCode = np.float
ctypedef np.float_t Real
# ctypedef np.ndarray[int, ndim=1] Index
# ctypedef vector[np.float_t] RealVector

cdef class LatticeSite
cdef class MacroBlock

# cdef class Index:
#     cdef vector[int]* v
#     def __cinit__(self):
#         self.v = new vector[int](3)

#     def __dealloc__(self):
#         del self.v

#     cdef inline object GetObject(Index self, np.ndarray array):
#         return array[deref(self.v)[0],
#                      deref(self.v)[1],
#                      deref(self.v)[2]]
#     cdef inline object SetObject(Index self, np.ndarray array, object val):
#         array[deref(self.v)[0],
#               deref(self.v)[1],
#               deref(self.v)[2]] = val
#         return val
    
#     cdef Real GetReal(Index self, np.ndarray array):
#         return array[deref(self.v)[0],
#                      deref(self.v)[1],
#                      deref(self.v)[2]]
    
cdef class Domain:
    """Represent the entire simulation domain that is needed for the
    simulation.
    """
    cdef:
        public Real VoxelSize
        public int BlockSize
        vector[Real]* _Origin
        vector[int]* _BlockCounts
        public np.ndarray blocks
        
    property Origin:
        def __get__(self):
            cdef np.ndarray ans = np.zeros(3, dtype=RealCode)
            for i from 0 <= i < 3:
                ans[i] = deref(self._Origin)[i]
            return ans
        def __set__(self, np.ndarray val):
            for i from 0 <= i < 3:
                deref(self._Origin)[i] = val[i]
    property BlockCounts:
        def __get__(self):
            cdef np.ndarray ans = np.zeros(3, dtype=int)
            for i from 0 <= i < 3:
                ans[i] = deref(self._BlockCounts)[i]
            return ans
        def __set__(self, np.ndarray val):
            for i from 0 <= i < 3:
                deref(self._BlockCounts)[i] = val[i]

    def __cinit__(self, *args, **kwargs):
        self._Origin = new vector[Real](3)
        self._BlockCounts = new vector[int](3)

    def __dealloc__(self):
        del self._Origin
        del self._BlockCounts
        
    def __init__(self, VoxelSize, SurfaceBounds, BlockSize=8):
        """VoxelSize - voxel size, in metres
        
        SurfaceBounds - bounds of the surface, in standard VTK order
        (x_min, x_max, y_min, y_max, z_min, z_max), in metres.
        """
        self.VoxelSize = VoxelSize
        self.BlockSize = BlockSize
        # VTK standard order of (x_min, x_max, y_min, y_max, z_min, z_max)
        bb = np.array(SurfaceBounds)
        bb.shape = (3, 2)

        origin = []
        blocks = []
        i = 0
        for min, max in bb:
            size = max - min
            # int() truncates, we add 2 to make sure there's enough
            # room for the sites just outside.
            nSites = int(size / VoxelSize) + 2

            # The extra space
            extra = nSites * VoxelSize - size
            # We want to balance this equally with the placement of
            # the first site.
            siteZero = min - 0.5 * extra

            nBlocks = nSites / BlockSize
            remainder = nSites % BlockSize
            if remainder:
                nBlocks += 1
                pass
            
            origin.append(siteZero)
            blocks.append(nBlocks)
            continue
        
        self.Origin = np.array(origin)
        self.BlockCounts = np.array(blocks)
        
        # Fill the blocks with Nones
        self.blocks = np.empty(self.BlockCounts, dtype=object)
        return

    cdef void CalcPositionFromIndex(self, vector[int]* index, vector[Real]* ans):
        for i from 0 <= i < 3:
            deref(ans)[i] = deref(self._Origin)[i] + self.VoxelSize * deref(index)[i]

    cdef MacroBlock GetBlock(self, np.ndarray[int] blockIjk):
        val = self.blocks[blockIjk[0], blockIjk[1], blockIjk[2]]
        if val is None:
            val = self.blocks[blockIjk[0], blockIjk[1], blockIjk[2]] = MacroBlock(self, blockIjk, self.BlockSize)
            pass
        return val

    cdef LatticeSite GetSite(self, np.ndarray[int] globalSiteIjk):
        cdef np.ndarray[int] blockIjk = globalSiteIjk / self.BlockSize
        cdef MacroBlock block = self.GetBlock(blockIjk)
        return block.GetSite(globalSiteIjk)

    pass

cdef class DomainSmartBlockIterator:
    cdef Domain domain
    cdef np.ndarray ijk, maxInds
    
    def __init__(self, Domain domain):
        self.domain = domain
        self.ijk = np.zeros(3, dtype=int)
        self.maxInds = domain.BlockCounts - 1
        self.ijk[2] = -1

    def __iter__(self): return self

    def __next__(self):
        cdef DomainSmartBlockIterator slf = self
        cdef np.ndarray[int] ijk = slf.ijk
        cdef int p, i, j, k
        cdef MacroBlock val
        
        # Delete any unnecessary blocks
        for i in range(ijk[0] - 1, ijk[0] + 1):
            if i < 0: continue
            if i == ijk[0] and i != slf.maxInds[0]: continue
            for j in range(ijk[1] - 1, ijk[1] + 1):
                if j < 0: continue
                if j == ijk[1] and j != slf.maxInds[1]: continue
                for k in range(ijk[2] - 1, ijk[2] + 1):
                    if k < 0: continue
                    if k == ijk[2] and k != slf.maxInds[2]: continue
                    slf.domain.blocks[i, j, k] = None
                    continue
                continue
            continue

        # Update the index vector
        ijk[2] += 1
        if ijk[2] == slf.domain.BlockCounts[2]:
            ijk[2] = 0
            
            ijk[1] += 1
            if ijk[1] == slf.domain.BlockCounts[1]:
                ijk[1] = 0

                ijk[0] += 1
                if ijk[0] == slf.domain.BlockCounts[0]:
                    raise StopIteration
        
        # If the block hasn't been created, do so
        val = slf.domain.blocks[ijk[0], ijk[1], ijk[2]]
        if val is None:
            val = slf.domain.blocks[ijk[0], ijk[1], ijk[2]] = MacroBlock(slf.domain, ijk, slf.domain.BlockSize)
        
        # "yield" the value
        return val

cdef class MacroBlock:
    cdef public Domain domain
    cdef public unsigned int size
    cdef vector[int]* ijk
    cdef public np.ndarray sites
    def __cinit__(self,  Domain domain, np.ndarray[int] ijk, unsigned int size):
        self.ijk = new vector[int](3)
        cdef int i
        for i in range(3):
            deref(self.ijk)[i] = ijk[i]

    def __dealloc__(self):
        del self.ijk
    
    def __init__(self, Domain domain, np.ndarray[int] ijk, unsigned int size):
        self.domain = domain
        self.size = size
        self.sites = np.empty((size, size, size), dtype=object)
        
        cdef np.ndarray[int] globalSiteIjk = np.zeros(3, dtype=int)
        # cdef vector[int] i = vector[int](3)
        cdef int i,j,k, dim
        
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    globalSiteIjk[0] = deref(self.ijk)[0] * size + i
                    globalSiteIjk[1] = deref(self.ijk)[1] * size + j
                    globalSiteIjk[2] = deref(self.ijk)[2] * size + k
                    self.sites[i,j,k] = LatticeSite(self, globalSiteIjk)
        return

    cpdef IterSites(self):
        return self.sites.flat

    cpdef NdEnumerateSites(self):
        return np.ndenumerate(self.sites)
    
    cdef LatticeSite GetLocalSite(self, np.ndarray[int] ijk):
        return self.sites[ijk[0], ijk[1], ijk[2]]

    cdef LatticeSite GetSite(MacroBlock self, np.ndarray[int] gIjk):
        cdef np.ndarray[int] lIjk = gIjk.copy()
        cdef int size = self.size
        cdef int dim
        for dim in range(3):
            lIjk[dim] -= deref(self.ijk)[dim]*size
        # Check if the coords belong to another block, i.e. any of
        # the local ones outside the range [0, self.size)
        if (lIjk[0]<0) or (lIjk[1]<0) or (lIjk[2]<0) or (lIjk[0]>=size) or (lIjk[1]>=size) or (lIjk[2]>=size):
            return self.domain.GetSite(gIjk)

        return self.GetLocalSite(lIjk)
    
    cdef void CalcPositionFromIndex(self, vector[int]* index, vector[Real]* ans):
        self.domain.CalcPositionFromIndex(index, ans)

cdef class LatticeSite:
    cdef public MacroBlock block
    cdef vector[int]* ijk
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
    
    def __cinit__(self, MacroBlock block, np.ndarray[int] ijk):
        self._Position = new vector[Real](3)
        self.block = block
        self.ijk = new vector[int](3)
        cdef int i
        for i in range(3):
            deref(self.ijk)[i] = ijk[i]
        
        block.domain.CalcPositionFromIndex(self.ijk, self._Position)

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
        
    def __dealloc__(self):
        del self._Position
        del self.ijk
        
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
    cdef Domain domain
    cdef np.ndarray shape
    cdef np.int_t i
    cdef np.int_t max_i
    cdef LatticeSite site
    cdef np.ndarray ind
    
    def __cinit__(self, LatticeSite site):
        self.site = site
        self.domain = site.block.domain
        self.shape = self.domain.BlockCounts * self.domain.BlockSize
        self.i = -1
        self.ind = np.zeros(3, dtype=int)
        return
    
    cdef np.ndarray GetVector(SiteNeighbourEnumeratorBase self):
        return neighbours[self.i]
    
    def __iter__(self): return self
    
    def __next__(self):
        cdef np.ndarray[int, ndim=1] latticeVec
        cdef SiteNeighbourEnumeratorBase slf = self
        cdef bint shouldTryNextVec = False
        while True:
            slf.i += 1
            if slf.i >= slf.max_i:
                raise StopIteration
            
            latticeVec = slf.GetVector()
            for j from 0 <= j < 3:
                slf.ind[j] = deref(slf.site.ijk)[j] + latticeVec[j]
                # Skip any out of range
                if slf.ind[j] < 0 or slf.ind[j] >= slf.shape[j]:
                    shouldTryNextVec = True
                    break
                shouldTryNextVec = False
                
            if shouldTryNextVec:
                continue
            else:
                break
        
        return slf.i, slf.site.block.GetSite(slf.ind)
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

