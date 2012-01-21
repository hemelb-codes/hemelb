#!/usr/bin/env python
import os.path
import numpy as np
import threading
import Queue
import collections

from hemeTools.parsers.geometry.generic import OutOfDomainSite, AllSolidBlock, Site, Block
from hemeTools.parsers.geometry.multiprocess import AsyncBlockProcessingLoader
from hemeTools.parsers import D3Q27Directions

import pdb

REALISTIC_IOLET_COUNT = 100

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

        errLines = [] 
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
        self.block = site.GetBlock()
        return
    pass

class CheckingLoader(AsyncBlockProcessingLoader):
    BLOCK_REPORT_PERIOD = 100
    
    def __init__(self, *args, **kwargs):
        self.Verbose = kwargs.pop('verbose', False)
        AsyncBlockProcessingLoader.__init__(self, *args, **kwargs)
        return

    def Info(self, message):
        if self.Verbose:
            print message
        return
    
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
            self.Info('Checked block {} of {}'
                      .format(bIjk, self.Domain.TotalBlocks))
            pass
        return AsyncBlockProcessingLoader.OnBlockProcessed(self,
                                                           bIdx)
        
    def OnEndPreamble(self):
        """Loaded the very basic domain information. Give a status
        report.
        """
        fields = ['Version', 'BlockSize', 'BlockCounts',
                  'VoxelSize', 'Origin']
        template = '\n'.join('%s: {0.%s}' % (f, f) for f in fields)
        self.Info(template.format(self.Domain))
        return
    
    def OnEndHeader(self):
        """All header data now loaded. Give a status report and check
        consistency of header data with itself and the actual file
        size.
        """
        fluidSiteCount = np.sum(self.Domain.BlockFluidSiteCounts)
        self.Info('NumberOfFluidSites: {}'.format(fluidSiteCount))
        
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
        hadErrors = False
        while True:
            # ReportQueue's a Queue.Queue so this is OK
            # (It will block until we get something)
            errList = self.ReportQueue.get()
            if errList is None:
                # Got the sentinal, so we need to exit
                # First, set a flag to indicate the presence of errors.
                
                # This assignment is OK since the main thread waits on
                # us after putting the sentinal into the queue (see
                # EndReportingThread).
                self.HadErrors = hadErrors
                break
            hadErrors = True
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
    
    def __call__(self, domain, bInd):
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
        block = domain.GetBlock(bInd)
        self.block = block
        blockErrors = BlockErrorCollection(block)
        self.numFluid = 0

        dom = block.GetDomain()
        
        if isinstance(block, AllSolidBlock):
            self.CheckFluidSiteCount(blockErrors.AddBlockError)
            return blockErrors.Format()
        
        # gs = global site
        # ls = local site
        for lsInd in np.ndindex((dom.BlockSize, dom.BlockSize, dom.BlockSize)):
            lsInd = np.array(lsInd)
            
            site = block.GetLocalSite(lsInd)

            # Quick helper function to add an error to this site's
            # error list
            addSiteError = lambda msg: blockErrors.AddSiteError(site, msg)
            
            self.CheckBasicFluidProperty(site, addSiteError)
            
            if not site.IsFluid:
                # Solid sites we don't check in detail since they have
                # no more data
                continue
            
            self.numFluid += 1
           
            self.CheckFluidSiteLinks(site, addSiteError)
            self.CheckIntersectionConsistency(site, addSiteError)
            self.CheckIOletConsistency(site, addSiteError)

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
    
    def CheckBasicFluidProperty(self, site, addSiteError):
        """Some basic sanity checking on the value of the 'isFluid' uint.
        """
        if not (site.IsFluid == Site.SOLID or
                site.IsFluid == Site.FLUID):
            addSiteError('Site doesn\'t appear to have an appropriate value for \'IsFluid\': %d'.format(site.IsFluid))
        return

    def CheckFluidSiteLinks(self, site, addSiteError):
        """Loop over the LB velocity set, checking link properties and
        whether the type (fluid only or fluid and edge) is consistent
        with the neighbours' types.  More concretely:

        - links to solid sites must have the intersection type set appropriately

        - links to fluid sites must have a NO_INTERSECTION, and vice versa
        """
        for iNeigh, delta in enumerate(D3Q27Directions):
            neigh = site.GetBlock().GetSite(site.Index + delta)

            if isinstance(neigh, OutOfDomainSite) or neigh.IsSolid:
                if site.IntersectionType[iNeigh] == Site.NO_INTERSECTION:
                    addSiteError('Fluid site has no intersection but neighbour is solid')
                pass
            else: # neigh is fluid
                if site.IntersectionType[iNeigh] != Site.NO_INTERSECTION:
                    addSiteError('Link to fluid site has the wrong kind of intersection (%d)'.format(site.IntersectionType[iNeigh]))
                
        return

    def CheckIntersectionConsistency(self, site, addSiteError):
        """Check that all links have a valid intersection type and 
        - that all links that aren't to other fluid sites have an intersection distance
        - that the distance (in lattice units) is in the range allowed
        """        
        for i in range(len(D3Q27Directions)):
            if site.IntersectionType[i] not in [Site.NO_INTERSECTION, Site.WALL_INTERSECTION, Site.INLET_INTERSECTION, Site.OUTLET_INTERSECTION]:
                addSiteError('Site had an invalid intersection type: %d'.format(site.IntersectionType[i]))
            if site.IntersectionType[i] != Site.NO_INTERSECTION:
                if site.IntersectionDistance[i] == None:
                    addSiteError('Site had a null intersection distance when it did have an intersection.')
                distance = np.sum(site.IntersectionDistance[i] ** 2) ** 0.5
                maxDistance = np.sum(D3Q27Directions[i] ** 2) ** 0.5
                if distance == 0 or distance >= maxDistance:
                    addSiteError('Site had an intersection distance (%f) outside the allowed range (0,%f)' % (distance, maxDistance))
        return
    
    def CheckIOletConsistency(self, site, addSiteError):
        """Walls are the solid boundaries of the simulation. Check
        some properties:
        
        - that the magnitude of the normal is very close to 1, and
        
        - that the distance (in lattice units) is in the range allowed
          for the lattice
        """
        for i in range(len(D3Q27Directions)):
            if site.IntersectionType[i] in [Site.INLET_INTERSECTION, Site.OUTLET_INTERSECTION]:
                if site.IOletIndex[i] < 0 or site.IOletIndex[i] > REALISTIC_IOLET_COUNT:
                    addSiteError('IOlet index for an iolet link (%d) was outside expected range.'.format(site.IOletIndex[i]))
            elif site.IOletIndex[i] != Site.NO_IOLET:
                addSiteError('Site had no link to an iolet site but DID have an iolet index set (%d).'.format(site.IOletIndex[i]))
        return
    pass
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Analyze a HemeLB input file for self-consistency.',
        epilog='The return value indicates the presence of detected errors.'
        )
    parser.add_argument('-q', '--quiet',
                        action='store_true', default=False,
                        help='enable quiet mode, which supresses non-error output')
    parser.add_argument('--debug',
                        action='store_true', default=False,
                        help='enable debugging output')
    parser.add_argument('input', nargs=1,
                        help='the input file to check')
    args = parser.parse_args()
    
    if args.debug:
        from hemeTools.parsers.config.multiprocess import GetLogger
        GetLogger().setLevel('DEBUG')
        pass
    
    ldr = CheckingLoader(args.input[0],
                         verbose=(not args.quiet))
    dom = ldr.Load()
    
    if ldr.HadErrors:
        raise SystemExit(1)
    
    


