# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import numpy as np
import xdrlib
from ...utils import xdr
import zlib
import os.path
from xml.etree import ElementTree

from .generic import Domain, Block, AllSolidBlock, Site
from .. import HemeLbMagicNumber
from . import GeometryMagicNumber

class GeometryParsingError(Exception):
    """Indicate an error while parsing a geometry file.
    """
    pass

class FakeUnpacker(object):
    """Fake xdrlib.Unpacker, for use when all the Sites in a Block are
    solid. Just returns zeros.
    """
    def unpack_uint(self):
        return Site.SOLID
    pass

PADDING_BYTE = 0

class ConfigLoader(object):
    """Loads a HemeLB config file.

    Offers a large number of hooks, in a manner similar to an
    event-driver XML parser for examining data on the fly.
    """
    VERSION = 4
    MIN_VERSION = 4
    MAX_VERSION = 4
    
    def __init__(self, filename):
        self.XmlFileName = filename
        tree = ElementTree.ElementTree()
        root = tree.parse(filename)
        # Get the GMY file from the XML
        gmyEl = root.find('geometry/datafile')
        assert gmyEl is not None, "XML does not have element 'geometry/datafile'!"
        relpath = gmyEl.get('path')
        assert relpath is not None, "Element 'geometry/datafile' does not have attribute 'path'"
        self.GmyFileName = os.path.join(os.path.dirname(filename), relpath)
        
        self.Domain = Domain()
        
        # Now the voxel size and origin
        vsEl = root.find('simulation/voxel_size')
        assert vsEl.get('units') == 'm', "voxel_size units not 'm'"
        self.Domain.VoxelSize = float(vsEl.get('value'))
        
        oEl = root.find('simulation/origin')
        assert oEl.get('units') == 'm', "origin units not 'm'"
        oStr = oEl.get('value')
        assert oStr[0] == '(' and oStr[-1] == ')'
        self.Domain.Origin = np.array(oStr[1:-1].split(','), dtype=float)
        assert self.Domain.Origin.shape == (3,)
        
        self.File = file(self.GmyFileName)
        return
    
    def Load(self):
        """Load the Domain encoded by the config file.
        """
        self._LoadPreamble()
        self._LoadHeader()
        self._LoadBody()
        return
    
    def _LoadPreamble(self):
        """Deal with the preamble (contains global stuff about the
        lattice: size, shape etc). Triggers OnBeginPreamble and
        OnEndPreamble events.
        1 uint for HemeLB magic number (0x686c6221)
        1 uint for geometry magic number (0x676d7904)
        1 uints for version number
        3 uints for domain size in blocks
        1 uint for number of sites along one side of a block
        1 uint of value 0 to pad to 72 bytes
        """
        self.OnBeginPreamble()
        
        self.PreambleBytes = 32
        preambleLoader = xdrlib.Unpacker(self.File.read(self.PreambleBytes))

        hlbNumber = preambleLoader.unpack_uint()
        if hlbNumber != HemeLbMagicNumber:
            raise GeometryParsingError("This doesn't appear to be a HemeLB file. Instead of the HemeLB magic number (%d), I found %d" % (HemeLbMagicNumber, hlbNumber))

        gmyNumber = preambleLoader.unpack_uint()
        if gmyNumber != GeometryMagicNumber:
            raise GeometryParsingError("This doesn't appear to be a geometry file. Instead of the geometry magic number (%d), I found %d" % (GeometryMagicNumber, gmyNumber))

        self.Domain.Version = preambleLoader.unpack_uint()
        if self.Domain.Version < self.MIN_VERSION:
            raise GeometryParsingError("This geometry file is of an older version (v%d) than this parser (v%d) can process" % (self.Domain.Version, self.VERSION))
        elif self.Domain.Version > self.MAX_VERSION:
            raise GeometryParsingError("This geometry file is of an newer version (v%d) than this parser (v%d) can process" % (self.Domain.Version, self.VERSION))
            pass
        
        self.Domain.BlockCounts = np.array([preambleLoader.unpack_uint() for i in xrange(3)], dtype=np.uint)
        self.Domain.BlockSize = preambleLoader.unpack_uint()
        
        padding = preambleLoader.unpack_uint()
        if padding != PADDING_BYTE:
            raise GeometryParsingError(r"The preamble to this file is padded with the wrong value. Instead of %d, I found %d" % (PADDING_BYTE, padding))

        self.OnEndPreamble()
        return

    def _LoadHeader(self):
        """Deal with the header containing per-block metadata (number
        of fluid sites and bytes occupied). Triggers OnBeginHeader and
        OnEndHeader events.
        """
        self.OnBeginHeader()
        self.HeaderBytes = self.Domain.TotalBlocks * 3 * 4

        headerLoader = xdrlib.Unpacker(self.File.read(self.HeaderBytes))
        
        nBlocks = self.Domain.TotalBlocks
        BlockCounts = self.Domain.BlockCounts
        
        BlockFluidSiteCounts = np.zeros(nBlocks, dtype=np.uint)
        BlockDataLength = np.zeros(nBlocks, dtype=np.uint)
        BlockUncompressedDataLength = np.zeros(nBlocks, dtype=np.uint)
        BlockStarts = np.zeros(nBlocks, dtype=np.uint)
        BlockEnds = np.zeros(nBlocks, dtype=np.uint)
        
        BlockEnds[-1] = self.PreambleBytes + self.HeaderBytes
        for bIjk in self.Domain.BlockIndexer.IterOne():
            BlockFluidSiteCounts[bIjk] = headerLoader.unpack_uint()
            BlockDataLength[bIjk] = headerLoader.unpack_uint()
            BlockUncompressedDataLength[bIjk] = headerLoader.unpack_uint()
            BlockStarts[bIjk] = BlockEnds[bIjk-1]
            BlockEnds[bIjk] = BlockStarts[bIjk] + BlockDataLength[bIjk]
            continue
        
        self.Domain.BlockFluidSiteCounts = BlockFluidSiteCounts
        self.BlockDataLength = BlockDataLength
        self.BlockUncompressedDataLength = BlockUncompressedDataLength
        self.BlockStarts = BlockStarts
        self.BlockEnds = BlockEnds

        self.OnEndHeader()
        return
    
    def _LoadBody(self):
        """Reads all the Blocks from the file, iterating in a z-index
        fastest fashion. Triggers OnBeginBody and OnEndBody.
        """
        self.OnBeginBody()
        
        dom = self.Domain
        
        nBlocks = dom.TotalBlocks
        BlockCounts = dom.BlockCounts
        
        for bIjk, bIdx in self.Domain.BlockIndexer.IterBoth():
            self._LoadBlock(dom, bIdx, bIjk)
            continue
        
        self.OnEndBody()
        
    def _LoadBlock(self, domain, bIdx, bIjk):
        """Loads a single Block. Iterating over Sites in a z-index
        fastest fashion.
        """
        self.OnBeginBlock(bIdx, bIjk)
        if domain.BlockFluidSiteCounts[bIjk] == 0:
            b = AllSolidBlock(domain, bIdx)
        else:
            b = Block(domain, bIdx)
        
            b.nFluidSites = domain.BlockFluidSiteCounts[bIjk]
            b.Sites = np.zeros(domain.BlockSize**3, dtype=object)
            
            compressed = self.File.read(self.BlockDataLength[bIjk])
            uncompressed = zlib.decompress(compressed)
            blockLoader = xdr.Unpacker(uncompressed)
            
            # sl = site local
            for slIjk, slIdx in domain.BlockSiteIndexer.IterBoth():
                # sg = site global
                sgIdx = domain.BlockSize * bIdx + slIdx
                b.Sites[slIjk] = self._LoadSite(b, sgIdx, blockLoader)
                continue
            pass
        
        domain.SetBlock(bIdx, b)
        
        self.OnEndBlock(bIdx, bIjk)
        return

    def _LoadSite(self, block, sgIdx, loader):
        """Loads a single site.
        """
        self.OnBeginSite(block, sgIdx)
        
        s = Site(block, sgIdx)
        s.LoadFrom(loader)
        self.OnEndSite(block, s)
        return s
    
    # Override these methods in your subclass.
    def OnBeginPreamble(self):
        return
    def OnEndPreamble(self):
        return
    def OnBeginHeader(self):
        return
    def OnEndHeader(self):
        return
    def OnBeginBody(self):
        return
    def OnEndBody(self):
        return
    def OnBeginBlock(self, bIdx, bIjk):
        return
    def OnEndBlock(self, bIdx, bIjk):
        return
    def OnBeginSite(self, block, sgIdx):
        return
    def OnEndSite(self, block, site):
        return
    
    pass
