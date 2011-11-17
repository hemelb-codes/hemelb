#!/usr/bin/env python
import os.path
import numpy as np
import threading
import Queue
import collections

from hemeTools.parsers.config.generic import OutOfDomainSite, AllSolidBlock
from hemeTools.parsers.config.cfg import *
from hemeTools.parsers.config import cfg
from hemeTools.parsers.config.multiprocess import AsyncBlockProcessingLoader

from multiprocessing.util import debug
import pdb

class EmptyErrorCollectionFormatError(Exception):
    pass

class Error(object):
    """Base class to hold info about the error in a config file.
    """
    def __init__(self, message):
        self.message = message
        return
    
    _template = '{message}'
    
    def Format(self):
        return self._template.format(**self.__dict__)
    
    def Indent(self, text):
        return '    ' + text
    
    pass
    
class DomainError(Error):
    """Config file errors to do with domain global data.
    """
    pass

class ErrorCollection(Error):
    """A container for Errors.
    """
    def __init__(self, item):
        self.item = item
        self.itemErrors = []
        self.subItemErrors = collections.OrderedDict()
        return

    @property
    def TotalErrors(self):
        return len(self.itemErrors) + len(self.subItemErrors)
    
    def Format(self, indent=0):
        nErr = self.TotalErrors
        if nErr == 0:
            return []
        elif nErr == 1:
            self.errString = 'Error'
        else:
            self.errString = 'Errors'
            pass
        
        errLines = ['{errString} in {item}:'.format(**self.__dict__)]
        for iErr in self.itemErrors:
            errLines.append(self.Indent(iErr.Format()))
            continue

        for siErrColl in self.subItemErrors.itervalues():
            try:
                subLines = siErrColl.Format()
                errLines.extend(self.Indent(sl) for sl in subLines)
            except EmptyErrorCollectionFormatError:
                pass
            continue

        return errLines
    pass

class BlockErrorCollection(ErrorCollection):
    """Container for errors to do with a particular block.
    """
    def AddBlockError(self, message):
        self.itemErrors.append(BlockError(self.item, message))
        return
    
    def AddSiteError(self, site, message):
        try:
            siteErrorCollection = self.subItemErrors[site]
        except KeyError:
            self.subItemErrors[site] = SiteErrorCollection(site)
            siteErrorCollection = self.subItemErrors[site]
            pass
        
        siteErrorCollection.Add(message)
        return
    
    pass
        
class SiteErrorCollection(ErrorCollection):
    """Container for errors to do with a particular site.
    """
    def Add(self, message):
        self.itemErrors.append(SiteError(self.item, message))
        return
    
    pass
    
class BlockError(Error):
    """An error attached a block's data.
    """
    def __init__(self, block, message):
        Error.__init__(self, message)
        self.block = block
        return
    pass

class SiteError(Error):
    """An error in a site's data.
    """
    def __init__(self, site, message):
        Error.__init__(self, message)
        self.site = site
        self.block = site.Block
        return
    pass

class D3Q15Lattice(object):
    """Represent needed properties of the LB lattice. Velocities must
    be ordered in the same way as the rest of HemeLB.
    """
    
    # Neighbour deltas are the vectors which when added to a fluid
    # site position vector give the location of a putative neighbour.
    velocities = np.array([[ 0,  0,  0],
                           [ 1,  0,  0],
                           [-1,  0,  0],
                           [ 0,  1,  0],
                           [ 0, -1,  0],
                           [ 0,  0,  1],
                           [ 0,  0, -1],
                           [ 1,  1,  1],
                           [-1, -1, -1],
                           [ 1,  1, -1],
                           [-1, -1,  1],
                           [ 1, -1,  1],
                           [-1,  1, -1],
                           [ 1, -1, -1],
                           [-1,  1,  1]])
    norms = np.sqrt(np.sum(velocities**2, axis=-1))
    minNorm = norms.min()
    maxNorm = norms.max()
    
    @classmethod
    def IsInNormRange(cls, x):
        return x >= cls.minNorm and x <= cls.maxNorm
    
    neighs = velocities[1:]
    pass

Lattice = D3Q15Lattice

class CheckingLoader(AsyncBlockProcessingLoader):
    BLOCK_REPORT_PERIOD = 100
    def Load(self):
        self.StartReportingThread()
        # Use a try finally to ensure the error thread gets cleaned
        # up, even if we get an exception.
        try:
            ans = AsyncBlockProcessingLoader.Load(self)
        finally:
            self.EndReportingThread()
            
        return ans

    def OnBlockProcessed(self, bIdx):
        bIjk = self.Domain.BlockIndexer.NdToOne(bIdx)
        if bIjk % self.BLOCK_REPORT_PERIOD == 0:
            print 'Checked block %d of %d' % (bIjk, self.Domain.TotalBlocks)
            pass
        return AsyncBlockProcessingLoader.OnBlockProcessed(self, bIdx)
        
    def OnEndPreamble(self):
        """Loaded the very basic domain information. Give a status
        report.
        """
        print ('-----\nInfo: stress type = ' +
               str(self.Domain.StressType) + ', sites per block = ' +
               str(self.Domain.BlockSize ** 3) + ', blocks: ' +
               str(self.Domain.BlockCounts[0]) + 'x' +
               str(self.Domain.BlockCounts[1]) + 'x' +
               str(self.Domain.BlockCounts[2]) + ', vox size = ' +
               str(self.Domain.VoxelSize) + ', origin = ' +
               str(self.Domain.Origin) + '\n-----\n')
        
        return
    
    def OnEndHeader(self):
        """All header data now loaded. Give a status report and check
        consistency of header data with itself and the actual file
        size.
        """
        fluidSiteCount = np.sum(self.Domain.BlockFluidSiteCounts)
        print '-----Fluid Site Count = %d\n-----\n' % fluidSiteCount
        
        # For consistency, if BlockDataLength[i] == 0 then
        # BlockFluidSiteCounts[i] must also be zero, and vice versa
        for bIjk, bIdx in self.Domain.BlockIndexer.IterBoth():
            if (self.BlockDataLength[bIjk] == 0 and
                self.Domain.BlockFluidSiteCounts[bIjk] != 0):
                
                self.PrintError(
                    BlockError(
                        self.Domain.GetBlock(bIdx),
                        'Header states no data but specifies some '
                        'fluid sites'
                        ).Format()
                    )
                pass
            
            if (self.Domain.BlockFluidSiteCounts[bIjk] == 0 and
                self.BlockDataLength[bIjk] != 0):
                self.PrintError(
                    BlockError(
                        self.Domain.GetBlock(bIdx),
                        'Header states no fluid sites but specifies '
                        'some data').Format()
                    )
                pass
            continue
        
        # The length of the file must be equal to the value we
        # calculate from the headers.
        claimedFileSize = (
            self.PreambleBytes + 
            self.HeaderBytes +
            np.sum(self.BlockDataLength)
            )
        if claimedFileSize != os.path.getsize(self.FileName):
            self.PrintError(
                DomainError(
                    'File length does not match file metadata'
                    ).Format()
                )
            pass
        
        self.Checker = BlockChecker()
        self.SetBlockProcessor(self.Checker)
        return

    def HandleBlockProcessingResult(self, result):
        # Result is a list of errors; only bother sending it if there
        # are items
        if len(result):
            self.ReportQueue.put(result)
        return
        
    def PrintError(self, err):
        self.ReportQueue.put([err])
        return
    
    def StartReportingThread(self):
        """Start a thread to deal with writing our reports.
        """
        self.ReportQueue = Queue.Queue()
        self.ReportingThread = threading.Thread(
            target=self.ReportingThreadMain,
            name='ReportingThread'
            )
        self.ReportingThread.start()
        return
    
    def ReportingThreadMain(self):
        """This is the entry point for the reporting thread. Make sure
        we're thread-safe here!

        We will just take items off the queue and print them unless
        they're the sentinal value (None), in which case we return,
        killing the reporing thread.
        """
        while True:
            errList = self.ReportQueue.get()
            if errList is None:
                # Got the sentinal
                break
            for e in errList:
                print e
    
    def EndReportingThread(self):
        """Kill the reporting thread by sending the sentinal value and
        clean up.
        """
        self.ReportQueue.put(None)
        self.ReportingThread.join()
        return
    
    pass

class BlockChecker(object):
    """Callable object that is run in a separate process to check a
    single block for consistency.
    
    """
    
    def __call__(self, block):
        """Now we iterate over every site on the block, and check
        constraints as specified in the Check* methods.
        
        A list of error STRINGS is returned, i.e. the result of
        calling Format() on the accumulated Error object.
        """

        # So this may LOOK dodgy, using an instance attribute to store
        # the errors when we've only instantiated one of these and are
        # checking a load of blocks in parallel. HOWEVER, we are
        # pickling the unexecuted instance created above out to each
        # process to be run per-block. Once results have been sent
        # back, the executed object is destroyed.
        self.block = block
        blockErrors = BlockErrorCollection(block)
        self.numFluid = 0
        
        dom = block.Domain
        
        if isinstance(block, AllSolidBlock):
            self.CheckFluidSiteCount(blockErrors.AddBlockError)
            return blockErrors.Format()
        
        bInd = block.Index
        
        # gs = global site
        # ls = local site
        for lsInd in np.ndindex((dom.BlockSize, dom.BlockSize, dom.BlockSize)):
            lsInd = np.array(lsInd)
            
            site = block.GetLocalSite(lsInd)
            # Quick helper function to add an error to this site's
            # error list
            addSiteError = lambda msg: blockErrors.AddSiteError(site, msg)
            
            self.CheckBasicConfigConsistency(site, addSiteError)
            
            if site.Type == cfg.SOLID_TYPE:
                # Solid sites we don't check in detail since they have
                # no more data
                continue
            
            self.numFluid += 1
           
            self.CheckFluidSiteLinks(site, addSiteError)
            self.CheckSiteBoundaryData(site, addSiteError)
            self.CheckSiteWallData(site, addSiteError)
            
            continue
        
        self.CheckFluidSiteCount(blockErrors.AddBlockError)
        return blockErrors.Format()
    
    def CheckFluidSiteCount(self, addBlockError):
        """The number of non-solid sites on each block must equal that
        declared in the header.
        """
        if self.numFluid != self.block.nFluidSites:
            addBlockError('File fluid site counts not equal to number loaded')
            pass
        return
    
    def CheckBasicConfigConsistency(self, site, addSiteError):
        """Some basic sanity checking on the value of the 'cfg'
        uint32.
        """
        if not (site.Config == cfg.SOLID_TYPE or
                site.Config == cfg.FLUID_TYPE):
            isTypeKnown = False
            if (site.Type == cfg.INLET_TYPE or
                site.Type == cfg.OUTLET_TYPE):
                isTypeKnown = True
                pass

            if site.IsEdge:
                isTypeKnown = True
                pass

            if not isTypeKnown:
                addSiteError('Site doesn\'t appear to have any fitting type for config = 0b{:032b}'.format(site.Config))
                pass
            pass
        return

    def CheckFluidSiteLinks(self, site, addSiteError):
        """Loop over the LB velocity set, checking link properties and
        whether the type (fluid only or fluid and edge) is consistent
        with the neighbours' types.  More concretely:

        - fluid sites must have no out-of-domain neighbours;

        - non-edge fluid sites must have no solid neighbours;

        - inlet/outlet/edge fluid sites must have a cut distance array;

        - links to solid sites must have a cut distance in [0,1];

        - links to fluid sites must either have no cut distance array
          or have an infinite cut distance, and

        - inlet/outlet/edge fluid sites must have at least one solid
          neighbour.

        """
        nSolidNeighs = 0
        for iNeigh, delta in enumerate(Lattice.neighs):
            neigh = site.Block.GetSite(site.Index + delta)

            if isinstance(neigh, OutOfDomainSite):
                addSiteError('Fluid site has out-of-domain neighbour {}'.format(neigh))
                pass
            
            if neigh.IsSolid:
                nSolidNeighs += 1

                if site.Type == cfg.FLUID_TYPE and not site.IsEdge:
                    addSiteError('Simple-fluid site has solid neighbour {}'.format(neigh))
                    pass

                if site.CutDistances is None:
                    addSiteError('No CutDistance array for fluid site with solid neighbour {}'.format(neigh))
                    
                elif (site.CutDistances[iNeigh] < 0. or
                    site.CutDistances[iNeigh] >= 1.):
                    
                    addSiteError('Link to solid {neigh} has invalid '
                                 'CutDistance = {cd}'.format(neigh=neigh,
                                                             cd=site.CutDistances[iNeigh]))
                    pass
                
            else:                       # neigh is fluid

                if (site.CutDistances is not None and 
                    site.CutDistances[iNeigh] != np.inf):
                    
                    addSiteError('Link to fluid {neigh} has non-infinite '
                                 'CutDistance = {cd}'.format(neigh=neigh,
                                                             cd=site.CutDistances[iNeigh]))
                    pass
                
                pass 

            continue                    # for loop over neighbours
        
        
        if ((site.Type == cfg.INLET_TYPE or
             site.Type == cfg.OUTLET_TYPE or
             site.IsEdge) and
            nSolidNeighs == 0):
            addSiteError('Non--simple-fluid site has all valid fluid neighbours.')
            pass
        
        return

    def CheckSiteBoundaryData(self, site, addSiteError):
        """Boundaries are inlets and outlets. Check some properties:
        
        - that the magnitude of the normal is very close to 1, and

        - that the distance is in [0, sqrt(3)]
        """
        
        if not (site.Type == cfg.INLET_TYPE or
                site.Type == cfg.OUTLET_TYPE):
            return
        
        normMag = np.sum(site.BoundaryNormal ** 2) ** 0.5
        
        if abs(normMag - 1) >= 0.001:
            addSiteError('Boundary normal {} has non-unity magnitude'.format(site.BoundaryNormal))
            pass
        
        if not Lattice.IsInNormRange(site.BoundaryDistance):
            addSiteError(
                'Boundary distance ({}) not in [{}, {}]'.format(
                    site.BoundaryDistance,
                    Lattice.minNorm,
                    Lattice.maxNorm
                    )
                )
            pass
        
        return
    
    def CheckSiteWallData(self, site, addSiteError):
        """Walls are the solid boundaries of the simulation. Check
        some properties:
        
        - that the magnitude of the normal is very close to 1, and

        - that the distance is in [0, sqrt(3)]
        """
        if not site.IsEdge:
            return
        
        normMag = np.sum(site.WallNormal ** 2) ** 0.5
        
        if abs(normMag - 1) >= 0.001:
            addSiteError('Wall normal {} has non-unity magnitude'.format(site.WallNormal))
            pass
        
        if not Lattice.IsInNormRange(site.WallDistance):
            addSiteError(
                'Wall distance ({}) not in [{}, {}]'.format(
                    site.WallDistance,
                    Lattice.minNorm,
                    Lattice.maxNorm
                    )
                )
            pass
        
        return
    pass
    
if __name__ == "__main__":
    import sys
    # Uncomment 3 lines below for debugging
    # import multiprocessing.util
    # logger = multiprocessing.log_to_stderr()
    # logger.setLevel(multiprocessing.SUBDEBUG)

    ldr = CheckingLoader(sys.argv[1])
    dom = ldr.Load()


    


