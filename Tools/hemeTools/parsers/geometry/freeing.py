# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import threading
import numpy as np

from .simple import ConfigLoader
from .generic import Block

class FreeingConfigLoader(ConfigLoader):
    """This ConfigLoader will ensure a Block is kept until all Blocks
    in it's neighbourhood (see below) have been marked as finished
    with through the OnBlockProcessed method.It will then delete the
    unneeded Block. Your subclass must ensure that OnBlockProcessed is
    called with the Idx of the finished Block once you have processed
    it. This is to enable asynchronous analysis of Blocks.
    
    A Block's neighbourhood is defined similarly to the mathematical
    concept of a neighbours; viz. as the Block itself and all the
    adjacent Blocks. These are the possible offsets to the
    neighbourhood.
    """
    NeighbourhoodOffsets = np.mgrid[-1:2, -1:2, -1:2].reshape((3,27)).transpose()

    def _LoadHeader(self):
        ConfigLoader._LoadHeader(self)
        bc = self.Domain.BlockCounts
        nBlocks = self.Domain.TotalBlocks
        # Number of blocks in its neighbourhood.
        self.BlockNeighbourhoodSize = np.zeros(nBlocks, dtype=np.uint8)
        # Number of blocks in the neighbourhood that are available
        self.BlockNeighbourhoodAvailable = np.zeros(nBlocks, dtype=np.uint8)
        # Number of blocks in the neighbourhood that are done
        self.BlockNeighbourhoodDone = np.zeros(nBlocks, dtype=np.uint8)
        # Is the block itself done
        self.IsBlockDone = np.zeros(nBlocks, dtype=np.bool)
        # Lock to ensure only one thread at a time updates IsBlockDone
        # and BlockNeighbourhoodDone
        self.DoneLock = threading.RLock()

        # Compute the size of each block's neighbourhood
        for bIjk, bIdx in self.Domain.BlockIndexer.IterBoth():
            for i, delta in enumerate(self.NeighbourhoodOffsets):
                nIdx = bIdx + delta
                if np.any(nIdx < 0) or np.any(nIdx >= bc):
                    continue
                
                nIjk = self.Domain.BlockIndexer.NdToOne(nIdx)
                self.BlockNeighbourhoodSize[nIjk] += 1
                continue
            continue
        return
    
    def _LoadBlock(self, domain, bIdx, bIjk):
        block = ConfigLoader._LoadBlock(self, domain, bIdx, bIjk)
        
        for delta in self.NeighbourhoodOffsets:
            nIdx = bIdx + delta
            if np.any(nIdx < 0) or np.any(nIdx >= self.Domain.BlockCounts):
                continue
            
            nIjk = self.Domain.BlockIndexer.NdToOne(nIdx)
            self.BlockNeighbourhoodAvailable[nIjk] += 1
            if self.BlockNeighbourhoodAvailable[nIjk] == self.BlockNeighbourhoodSize[nIjk]:
                self.OnBlockNeighboursAvailable(nIdx)
                
            continue
        
        return

    def OnBlockNeighboursAvailable(self, bIdx):
        """Override this method to trigger processing of a block that
        may require its neighbours to be present.  When you have
        finished processing, you must call OnBlockProcessed with the
        3D index of the finished block.
        """
        return

    def OnBlockProcessed(self, bIdx):
        """Mark a block as finished with.
        """
        bIjk = self.Domain.BlockIndexer.NdToOne(bIdx)
        with self.DoneLock:
            self.IsBlockDone[bIjk] = True
            for delta in self.NeighbourhoodOffsets:
                nIdx = bIdx + delta
                if np.any(nIdx < 0) or np.any(nIdx >= self.Domain.BlockCounts):
                    continue

                nIjk = self.Domain.BlockIndexer.NdToOne(nIdx)
                self.BlockNeighbourhoodDone[nIjk] += 1
                if self.BlockNeighbourhoodDone[nIjk] == self.BlockNeighbourhoodSize[nIjk]:
                    self.OnBlockNeighboursProcessed(nIdx)
                    pass
                continue
            if np.alltrue(self.IsBlockDone):
                self.OnAllBlocksProcessed()
                
        return

    def OnAllBlocksProcessed(self):
        """This method is triggered once all blocks have been marked
        as processed.
        """
        return
    
    def OnBlockNeighboursProcessed(self, bIdx):
        self.Domain.DeleteBlock(bIdx)
        return
    
    pass
