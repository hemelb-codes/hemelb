# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import itertools
import numpy as np
import weakref

class NdIndexConverter(object):
    """Help for converting between 1d and Nd indices into arrays and
    iterating over them. This assumes that the arrays have a
    C-language ordering, i.e. last-fastest.
    """
    
    def __init__(self, shape):
        """shape is the number of elements along each dimension of the
        array.  """
        
        self.shape = np.array(shape).squeeze()
        self.ndim = len(self.shape)
        self._conv = np.array([np.prod(self.shape[i:]) for i in range(1, self.ndim+1)], dtype=np.int)
        return
    
    def OneToNd(self, one):
        """Go from a 1d index to an Nd index.
        """
        ans = np.zeros(self.ndim, dtype=np.int)
        for i in xrange(self.ndim):
            ans[i] = one / self._conv[i]
            one = one % self._conv[i]
            continue
        return ans

    def NdToOne(self, nd):
        """Go from an Nd index list to a 1d index.
        """
        return np.dot(nd, self._conv)

    def IterOne(self):
        """Return an iterator over the 1d indices.
        """
        return xrange(np.prod(self.shape))
    
    def IterNd(self):
        """Return an iterator over the Nd indices.
        """
        return itertools.imap(np.array, np.ndindex(tuple(self.shape)))

    def IterBoth(self):
        """Return an iterator yielding (1d, Nd) index tuples.
        """
        return itertools.izip(self.IterOne(), self.IterNd())
    
    pass

# ======================================================================
# Classes representing HemeLB objects
# ======================================================================

class Domain(object):
    """Spatial domain sampled by a HemeLB config file.

    Public attributes:

    - HemeLB magic number
    - Geometry magic number
    - Version number
    - BlockCounts
    - BlockSize
    - VoxelSize
    - Origin
    - SiteCounts
    - TotalBlocks

    Public methods:

    - GetBlock(bIdx)
    - SetBlock(bIdx, block)
    - GetSite(sgIdx)
    
    """
    def __init__(self):
        self._Version = None
        self._BlockCounts = None
        self._BlockSize = None
        
        self.VoxelSize = None
        self.Origin = None
        
        return

    @property
    def SiteCounts(self):
        return self.BlockCounts * self.BlockSize

    @property
    def BlockCounts(self):
        return self._BlockCounts
    
    @BlockCounts.setter
    def BlockCounts(self, bc):
        self._BlockCounts = bc
        self.Blocks = np.zeros(self.TotalBlocks, dtype=object)
        self.BlockIndexer = NdIndexConverter(bc)
        for ijk, i3 in self.BlockIndexer.IterBoth():
            self.Blocks[ijk] = NotYetLoadedBlock(self, i3)
            continue
        
        return
        
    @property
    def TotalBlocks(self):
        # Numpy up-casting rules mean that numpy.uint x python int
        # may promote to a float, but this is platform dependent...
        # Ensure that this is python int to be platform independent.
        return int(np.prod(self._BlockCounts))

    @property
    def BlockSize(self):
        return self._BlockSize
    @BlockSize.setter
    def BlockSize(self, bs):
        self._BlockSize = bs
        self.BlockSiteIndexer = NdIndexConverter((bs,bs,bs))
        return
    
    def GetBlock(self, bIdx):
        if np.all(bIdx >= 0) and np.all(bIdx < self.BlockCounts):
            return self.Blocks[self.BlockIndexer.NdToOne(bIdx)]
        else:
            return OutOfDomainBlock(self, bIdx)
        
    def SetBlock(self, bIdx, block):
        if np.all(bIdx >= 0) and np.all(bIdx < self.BlockCounts):
            self.Blocks[self.BlockIndexer.NdToOne(bIdx)] = block
        else:
            raise ValueError('Block index out of range')

    def DeleteBlock(self, bIdx):
        """Ensure that the circular references from Sites to their
        containing Block are deleted.
        """
        self.GetBlock(bIdx).DeleteSites()
        self.SetBlock(bIdx, None)

    def GetSite(self, sIdx):
        bIdx = sIdx / self.BlockSize
        block = self.GetBlock(bIdx)
        return block.GetLocalSite(sIdx % self.BlockSize)

    pass

class Block(object):
    """A macro-block, typically 8x8x8, but this is variable here.

    Public attributes:

     - Domain
     - Index

    Public methods:

     - GetSite(sgIdx)
    """
    def __init__(self, domain, index):
        self.GetDomain = weakref.ref(domain)
        self.Index = index
        self.nFluidSites = 0
        self.Sites = None
        return
    
    def __getstate__(self):
        picdic = self.__dict__.copy()
        picdic['Domain'] = self.GetDomain()
        assert picdic['Domain'] is not None
        del picdic['GetDomain']
        return picdic
    
    def __setstate__(self, picdic):
        picdic['GetDomain'] = weakref.ref(picdic['Domain'])
        del picdic['Domain']
        self.__dict__.update(picdic)
        return

    def GetLocalSite(self, slIdx):
        assert np.all(slIdx >= 0) and np.all(slIdx < self.GetDomain().BlockSize)
        return self.Sites[self.GetDomain().BlockSiteIndexer.NdToOne(slIdx)]

    def GetSite(self, sgIdx):
        return self.GetDomain().GetSite(sgIdx)

    def DeleteSites(self):
        del self.Sites
        return
    
    _template = 'Block [' + ', '.join('{0[%d]:{2[%d]}}/{1[%d]:{2[%d]}}' % (i,i,i,i) for i in xrange(3)) + ']'
    
    def __format__(self, format_spec):
        bc = self.GetDomain().BlockCounts
        # Add something less than one, to ensure integer powers of ten are 
        # rounded up.
        widths = np.ceil(np.log10(bc + 0.1)).astype(int)
        
        return self._template.format(self.Index, bc, widths)
    
    pass

    
class NotYetLoadedBlock(Block):
    def GetLocalSite(self, sIdx):
        raise ValueError('Cannot get sites from NotYetLoadedBlock')
    pass

class Site(object):
    DIRECTIONS = 26

    SOLID = 0
    FLUID = 1

    NO_IOLET = -1

    NO_INTERSECTION = 0
    WALL_INTERSECTION = 1
    INLET_INTERSECTION = 2
    OUTLET_INTERSECTION = 3
    INTERSECTION_TYPES = (NO_INTERSECTION,
                          WALL_INTERSECTION,
                          INLET_INTERSECTION,
                          OUTLET_INTERSECTION)
    def __init__(self, block, sgIdx):
        self.GetBlock = weakref.ref(block)
        self.Index = sgIdx
        self.Type = Site.SOLID
        
        self.IntersectionType = None
        self.IntersectionDistance = None
        self.IOletIndex = None

        dom = block.GetDomain()
        self.Position = dom.Origin + dom.VoxelSize * sgIdx
        
        return
    
    def __getstate__(self):
        picdic = self.__dict__.copy()
        picdic['Block'] = self.GetBlock()
        assert picdic['Block'] is not None
        del picdic['GetBlock']
        return picdic
    
    def __setstate__(self, picdic):
        picdic['GetBlock'] = weakref.ref(picdic['Block'])
        del picdic['Block']
        self.__dict__.update(picdic)
        return
    

    @property
    def IsEdge(self):
        return bool(self.IsFluid and 
                    np.any(self.IntersectionType == self.WALL_INTERSECTION))

    @property
    def IsFluid(self):
        return self.Type == Site.FLUID
    
    @property
    def IsSolid(self):
        return self.Type == Site.SOLID
    
    _template = 'Site [' + ', '.join('{0[%d]:{2[%d]}}/{1[%d]:{2[%d]}}' % (i,i,i,i) for i in xrange(3)) + ']'
    
    def __format__(self, format_spec):
        sc = self.GetBlock().GetDomain().SiteCounts
        widths = np.ceil(np.log10(sc)).astype(int)
        
        return self._template.format(self.Index, sc, widths)
    
    pass

class OutOfDomainBlock(Block):
    
    def GetLocalSite(self, slIdx):
        assert np.all(slIdx >= 0) and np.all(slIdx < self.GetDomain().BlockSize)
        sgIdx = self.Index * self.GetDomain().BlockSize + slIdx
        return OutOfDomainSite(self, sgIdx)
    
    pass

class OutOfDomainSite(Site):
    def __init__(self, block, sgIdx):
        self.IsFluid = False
        self.GetBlock = weakref.ref(block)
        self.Index = sgIdx
        return
    pass


class AllSolidBlock(Block):
    def GetLocalSite(self, slIndx):
        assert np.all(slIndx >= 0) and np.all(slIndx < self.GetDomain().BlockSize)
        return AllSolidSite(self, slIndx)
    pass

class AllSolidSite(Site):
    def __init__(self, block, sgIdx):
        self.GetBlock = weakref.ref(block)
        self.Index = sgIdx
        self.Type = Site.SOLID
        return
    pass

