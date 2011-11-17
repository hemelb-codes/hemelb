"""ConfigLoader which will process blocks asynchronously in parallel.

For debugging, add the following to the top of your __main__:

import multiprocessing
logger = multiprocessing.log_to_stderr()
logger.setLevel(multiprocessing.SUBDEBUG)

and use multiprocessing.util.debug to print messages within the code
to be executed per-block. This will also print tracebacks from any
uncaught exceptions from your per-block code.
"""
import multiprocessing
import copy
import signal
import traceback
import numpy as np
import pdb
from .generic import NdIndexConverter, Domain, Block, OutOfDomainBlock, AllSolidBlock
from .freeing import FreeingConfigLoader

from multiprocessing.util import debug

class SubDomain(Domain):
    """A Domain subclass that only has a subset of the full
    domain. This is useful to avoid sending very large messages over
    interprocess Queues.
    """
    
    def __init__(self, parent, bIdx):
        """Constructor; requires the domain it's mimicking and the 3D
        index of the block at the centre.
        """
        self.StressType = parent.StressType
        self._BlockCounts = parent._BlockCounts
        self._BlockSize = parent._BlockSize
        
        self.VoxelSize = parent.VoxelSize
        self.Origin = parent.Origin
        
        self.BlockIndexer = parent.BlockIndexer
        self.BlockSiteIndexer = parent.BlockSiteIndexer

        self._Index = bIdx
        self._Shift = bIdx - 1
        
        self.Blocks = np.empty(27, dtype=object)
        self.LocalBlockIndexer = NdIndexConverter((3,3,3))
        for lbIjk, lbIdx in self.LocalBlockIndexer.IterBoth():
            pbIdx = lbIdx + self._Shift
            pb = parent.GetBlock(pbIdx)
            if isinstance(pb, OutOfDomainBlock):
                lb = OutOfDomainBlock(self, pbIdx)
            elif isinstance(pb, AllSolidBlock):
                lb = AllSolidBlock(self, pbIdx)
            else:
                lb = PartialBlock(self, pb, lbIdx - 1)
            self.Blocks[lbIjk] = lb
            continue
        
        return
    
    def GetBlock(self, pbIndex):
        # The block index is from the parent Domain's array, must
        # convert to our local index.
        lbIndex = pbIndex - self._Shift
        if np.all(lbIndex >= 0) and np.all(lbIndex < 3):
            return self.Blocks[self.LocalBlockIndexer.NdToOne(lbIndex)]
        else:
            raise ValueError("Don't have that block")
        
    def SetBlock(self, bIndex, block):
        raise ValueError("Can't set blocks!")
        
    def GetSite(self, sIndex):
        bIndex = sIndex / self.BlockSize
        block = self.GetBlock(bIndex)
        return block.GetLocalSite(sIndex % self.BlockSize)

    pass

class PartialBlock(Block):
    """A Block that only actually has sites for the parts needed by
    the Block at the centre of the SubDomain.  Useful to further
    reduce data transferring over IPC. Viz. a block with a <100> type
    offset from the centre only has one face of Sites, <110> offset
    has only one edge and a <111> offset has only a single corner.
    """
    def __init__(self, newDomain, parent, targetToThis):
        """Constructor; requires the SubDomain this is part of, the
        Block this is mirroring and the vector from the block at the
        centre of the SubDomain to this Block.
        """
        self.Domain = newDomain
        self.Index = parent.Index
        self.nFluidSites = parent.nFluidSites
        
        self.Sites = np.empty(parent.Sites.shape, dtype=object)

        # We only want to copy the bits we care about
        
        # Here, set up iterators for the accessible subset of of this
        # block, from the target block
        bs = newDomain.BlockSize
        iters = []
        for diff in targetToThis:
            if diff < 0:
                it = xrange(bs - 1, bs)
            elif diff == 0:
                it = xrange(0, bs)
            else:
                it = xrange(0, 1)
                pass
            iters.append(it)
            continue
        
        # Now iterate over the interesting bit
        sIdx = np.zeros(3, dtype=np.int)
        for sIdx[0] in iters[0]:
            for sIdx[1] in iters[1]:
                for sIdx[2] in iters[2]:
                    sIjk = newDomain.BlockSiteIndexer.NdToOne(sIdx)
                    site = copy.copy(parent.Sites[sIjk])
                    site.Block = self
                    self.Sites[sIjk] = site
                    continue
                continue
            continue
        
        return

    def GetLocalSite(self, sIndex):
        assert np.all(sIndex >= 0) and np.all(sIndex < self.Domain.BlockSize)
        return self.Sites[self.Domain.BlockSiteIndexer.NdToOne(sIndex)]

    def GetSite(self, sIndex):
        return self.Domain.GetSite(sIndex)
    
    pass

class AsyncBlockProcessingLoader(FreeingConfigLoader):
    """A ConfigLoader that will asynchronously and in parallel apply a
    processing function to each Block in the Domain.

    You must make the following attributes callables with the
    specified signatures in your subclass.

    You must assign thProcessBlock(block) -> result
    
    HandleBlockProcessingResult(result) -> None

    This class will ensure that OnBlockProcessed is called
    appropriately once HandleBlockProcessingResult has returned.
    
    """
    @staticmethod
    def InitializeWorker():
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        return
    
    def __init__(self, *args, **kwargs):
        self._Workers = multiprocessing.Pool(initializer=self.InitializeWorker)
        FreeingConfigLoader.__init__(self, *args, **kwargs)
        return
    
    def SetBlockProcessor(self, bp):
        """Specify the object to process each Block. This MUST be a
        pickable object, unfortunately instance methods are not
        pickable. Further, this object will be serialised and sent to
        the worker processes so if should be as small as possible.
        """
        self._BlockProcessor = BlockProcessorWrapper(bp)
        return
    
    def OnBlockNeighboursAvailable(self, bIdx):
        debug('Submit block ' + str(bIdx))
        sd = SubDomain(self.Domain, bIdx)
        b = sd.GetBlock(bIdx)
        self._Workers.apply_async(self._BlockProcessor,
                                  args=(b, ),
                                  callback=self._HandleBlockResult)
        return
    
    def _HandleBlockResult(self, args):
        bIdx, result = args
        debug('Got results for block ' + str(bIdx))
        self.HandleBlockProcessingResult(result)
        self.OnBlockProcessed(bIdx)
        return

    def _LoadBody(self):
        # This try block ensures that the Pool is cleaned propery. In
        # the normal case, the input queue is closed and drained. If
        # an exception occurs, the workers are killed before the
        # exception is allow to propagate. In both cases the processes
        # must be joined.
        try:
            FreeingConfigLoader._LoadBody(self)
        except:
            self._Workers.terminate()
            raise
        else:
            self._Workers.close()
        finally:
            self._Workers.join()
        return

    pass

class BlockProcessorWrapper(object):
    """Wrapper around the Block-processing callable.

    This ensures the result is packaged up with the block index for
    AsyncBlockProcessingLoader._HandleBlockResults above and ensures
    any uncaught exceptions are logged.
    """
    
    def __init__(self, callable):
        self.callable = callable
        return
    
    def __call__(self, block):
        """Executed in the worker process. This deals with calling the
        user function and returning the result packaged up with the
        block index.
        """
        
        debug('Process block %s' % str(block.Index))
        try:
            result = self.callable(block)
            
        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            debug(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise
        
        # It was fine, give a normal answer
        return (block.Index, result)
    
    pass
