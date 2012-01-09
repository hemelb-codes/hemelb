import numpy as np
import xdrlib

from . import cfg
from .generic import Domain, Block, AllSolidBlock, Site

class FakeUnpacker(object):
    """Fake xdrlib.Unpacker, for use when all the Sites in a Block are
    solid. Just returns zeros.
    """
    
    def unpack_uint(self):
        return cfg.SOLID_TYPE
    pass

class ConfigLoader(object):
    """Loads a HemeLB config file.

    Offers a large number of hooks, in a manner similar to an
    event-driver XML parser for examining data on the fly.
    """
    
    def __init__(self, filename):
        self.FileName = filename
        self.File = file(filename)
        self.Domain = d = Domain()
        return
    
    def Load(self):
        """Load and return the Domain encoded by the config file.
        """
        self._LoadPreamble()
        self._LoadHeader()
        self._LoadBody()
        return self.Domain
    
    def _LoadPreamble(self):
        """Deal with the preamble (contains global stuff about the
        lattice: size, shape etc). Triggers OnBeginPreamble and
        OnEndPreamble events.
        """
        self.OnBeginPreamble()
        
        self.PreambleBytes = 4*4 + 4*8
        preambleLoader = xdrlib.Unpacker(self.File.read(self.PreambleBytes))
        
        self.Domain.BlockCounts = np.array([preambleLoader.unpack_uint() for i in xrange(3)], dtype=np.uint)
        self.Domain.BlockSize = preambleLoader.unpack_uint()
        
        self.Domain.VoxelSize = preambleLoader.unpack_double()
        self.Domain.Origin = np.array([preambleLoader.unpack_double() for i in xrange(3)], dtype=np.double)
        
        self.OnEndPreamble()
        return

    def _LoadHeader(self):
        """Deal with the header containing per-block metadata (number
        of fluid sites and bytes occupied). Triggers OnBeginHeader and
        OnEndHeader events.
        """
        self.OnBeginHeader()
        self.HeaderBytes = self.Domain.TotalBlocks * 2 * 4

        headerLoader = xdrlib.Unpacker(self.File.read(self.HeaderBytes))
        
        nBlocks = self.Domain.TotalBlocks
        BlockCounts = self.Domain.BlockCounts
        
        BlockFluidSiteCounts = np.zeros(nBlocks, dtype=np.uint)
        BlockDataLength = np.zeros(nBlocks, dtype=np.uint)
        BlockStarts = np.zeros(nBlocks, dtype=np.uint)
        BlockEnds = np.zeros(nBlocks, dtype=np.uint)
        
        BlockEnds[-1] = self.PreambleBytes + self.HeaderBytes
        for bIjk in self.Domain.BlockIndexer.IterOne():
            BlockFluidSiteCounts[bIjk] = headerLoader.unpack_uint()
            BlockDataLength[bIjk] = headerLoader.unpack_uint()
            BlockStarts[bIjk] = BlockEnds[bIjk-1]
            BlockEnds[bIjk] = BlockStarts[bIjk] + BlockDataLength[bIjk]
            continue
        
        self.Domain.BlockFluidSiteCounts = BlockFluidSiteCounts
        self.BlockDataLength = BlockDataLength
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
        
            blockLoader = xdrlib.Unpacker(self.File.read(self.BlockDataLength[bIjk]))
            
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
        
        s.Config = loader.unpack_uint()
        # Solid and simple fluid, we are done loading
        if not (s.Config == cfg.SOLID_TYPE or s.Config == cfg.FLUID_TYPE):
            
            if s.Type == cfg.INLET_TYPE or s.Type == cfg.OUTLET_TYPE:
                s.BoundaryNormal = np.array([loader.unpack_double() for i in xrange(3)], dtype=np.double)
                s.BoundaryDistance = loader.unpack_double()
                            
            if s.IsEdge:
                s.WallNormal = np.array([loader.unpack_double() for i in xrange(3)], dtype=np.double)
                s.WallDistance = loader.unpack_double()
                                        
            s.CutDistances = np.array([loader.unpack_double() for i in xrange(14)], dtype=np.double)
            pass
        
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
