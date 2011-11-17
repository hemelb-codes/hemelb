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

    neighs = velocities[1:]
    
    def OnEndPreamble(self):
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
        print '-----Fluid Site Count = %d\n-----\n' % np.sum(self.Domain.BlockFluidSiteCounts)
        # CONSTRAINT 0: if BlockDataLength[i] == 0 then BlockFluidSiteCounts[i] must also be zero and vice versa
        for bIjk, bIdx in self.Domain.BlockIndexer.IterBoth():
            if self.BlockDataLength[bIjk] == 0 and self.Domain.BlockFluidSiteCounts[bIjk] != 0:
                self.PrintError(
                    BlockError(self.Domain.GetBlock(bIdx),
                                'Header states no data but specifies some fluid sites').Format()
                                )
                
            if self.Domain.BlockFluidSiteCounts[bIjk] == 0 and self.BlockDataLength[bIjk] != 0:
                self.PrintError(
                    BlockError(self.Domain.GetBlock(bIdx),
                                'Header states no fluid sites but specifies some data').Format()
                                )
                
                    
        # CONSTRAINT 1: the length of the file must be equal to the
        # value we calculate from the headers.
        if (self.PreambleBytes + self.HeaderBytes + np.sum(self.BlockDataLength)) != os.path.getsize(self.FileName):
            self.PrintError(DomainError('File length does not match file metadata').Format())
            
        self.Checker = BlockChecker(self)
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
        self.ReportingThread = threading.Thread(target=self.ReportingThreadMain,
                                                    name='ReportingThread')
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
    def __init__(self, loader):        
        self.neighs = loader.neighs
        self.domainExtremeSite = loader.Domain.BlockCounts * loader.Domain.BlockSize - 1
        return

    def __call__(self, block):
        """Now we iterate over every site on the block, and check the following constraints.
        
        CONSTRAINT 2: the number of non-solid sites on each block must equal that declared in the header.
        CONSTRAINT 3: every site which has a fluid type (and isn't an edge) can have no invalid or solid-type neighbours.
        CONSTRAINT 4: every site which is not a fluid type must have at least one invalid or solid-type neighbour.
        CONSTRAINT 5: if the site is not of fluid or solid type, at least one of the cut distances must be between 0 and 1
        CONSTRAINT 6: if the site is not of fluid or solid type, all cut distances must be between 0 and 1 or be infinite
        CONSTRAINT 7: if the site is of inlet or outlet type, the magnitude of the boundary normal must be 1
        CONSTRAINT 8: if the site is of inlet or outlet type, the boundary distance must be between 0 and 1
        CONSTRAINT 9: if the site is an edge, the magnitude of the wall normal must be 1
        CONSTRAINT 10: if the site is an edge, the distance to the wall must be between 0 and root(3)

        A list of error messages is returned.
        """

        # So this may look dodgy, using an instance attribute to store
        # the errors when we are running a load of these in
        # parallel. HOWEVER, we are pickling the unexecuted instance
        # created above out to each process to be run serially.
        self.block = block
        self.blockErrors = BlockErrorCollection(block)
        self.numFluid = 0
        
        dom = block.Domain
        
        if isinstance(block, AllSolidBlock):
            self.CheckFluidSiteCount()
            return self.blockErrors.Format()
        
        bInd = block.Index
        
        # gs = global site
        # ls = local site
        for lsInd in np.ndindex((dom.BlockSize, dom.BlockSize, dom.BlockSize)):
            lsInd = np.array(lsInd)
            
            site = block.GetLocalSite(lsInd)
            # quick helper function to add an error to this site's error list
            AddSiteError = lambda msg: self.blockErrors.AddSiteError(site, msg)

            # Some basic sanity checking on the cfg
            if not (site.Config == cfg.SOLID_TYPE or site.Config == cfg.FLUID_TYPE):
                isTypeKnown = False
                if site.Type == cfg.INLET_TYPE or site.Type == cfg.OUTLET_TYPE:
                    isTypeKnown = True
                    pass

                if site.IsEdge:
                    isTypeKnown = True
                    pass

                if not isTypeKnown:
                    AddSiteError('Site doesn\'t appear to have any fitting type for config = 0b{:032b}')
                    pass
                pass

            if site.Type == cfg.SOLID_TYPE:
                # Solid sites we don't check in detail since they have no more data
               continue
            else:
                self.numFluid += 1
                pass

            # Loop over the LB velocity set, checking link properties
            nSolidNeighs = 0
            for iNeigh, delta in enumerate(self.neighs):
                neigh = block.GetSite(site.Index + delta)

                if isinstance(neigh, OutOfDomainSite):
                    AddSiteError('Fluid site has out-of-bounds neighbour {}'.format(neigh))
                if neigh.IsSolid:
                    nSolidNeighs += 1

                    if site.Type == cfg.FLUID_TYPE and not site.IsEdge:
                        AddSiteError('Simple-fluid site has solid neighbour {}'.format(neigh))
                        pass

                    if site.CutDistances[iNeigh] < 0. or site.CutDistances[iNeigh] >= 1.:
                        AddSiteError('Link to solid {neigh} has invalid CutDistance = {cd}'.format(neigh=neigh, cd=site.CutDistances[iNeigh]))

                else:           # neigh is fluid

                    if site.CutDistances is not None and \
                        site.CutDistances[iNeigh] != np.inf:
                        AddSiteError('Link to fluid {neigh} has non-infinite CutDistance = {cd}'.format(neigh=neigh, cd=site.CutDistances[iNeigh]))
                        pass

                    pass

                continue        # for loop over neighbours

            # CON 4
            if (site.Type == cfg.INLET_TYPE or site.Type == cfg.OUTLET_TYPE or site.IsEdge) and nSolidNeighs == 0:
                    AddSiteError('Non--simple-fluid site has all valid fluid neighbours.')
                    pass

            if site.Type == cfg.INLET_TYPE or site.Type == cfg.OUTLET_TYPE:
                normMag = np.sum(site.BoundaryNormal ** 2) ** 0.5
                # CON 7
                if abs(normMag - 1) >= 0.001:
                    AddSiteError('Boundary normal {} has non-unity magnitude'.format(site.BoundaryNormal))
                # CON 8
                if site.BoundaryDistance <= 0. or site.BoundaryDistance >= np.sqrt(3.):
                    AddSiteError('Boundary distance ({}) not in [0,sqrt(3)]'.format(site.BoundaryDistance))

            if site.IsEdge:
                normMag = np.sum(site.WallNormal ** 2) ** 0.5
                # CON 9
                if abs(normMag - 1) >= 0.001:
                    AddSiteError('Wall normal {} has non-unity magnitude'.format(site.WallNormal))
                # CON 10
                if site.WallDistance <= 0. or site.WallDistance >= (3.0**0.5):
                            AddSiteError('Wall distance ({}) not in [0,sqrt(3)]'.format(site.WallDistance))
            continue
        self.CheckFluidSiteCount()
        return self.blockErrors.Format()
    
    def CheckFluidSiteCount(self):
        # CON 2
        if self.numFluid != self.block.nFluidSites:
            self.blockErrors.AddBlockError('File fluid site counts not equal to number loaded')
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


    


