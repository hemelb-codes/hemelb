#!/usr/bin/env python
import os.path
import numpy as np
import threading
import Queue

from hemeTools.parsers.config.freeing import FreeingConfigLoader
from hemeTools.parsers.config.generic import AllSolidBlock
from hemeTools.parsers.config.cfg import *
from hemeTools.parsers.config import cfg
from hemeTools.parsers.config.multiprocess import AsyncBlockProcessingLoader
import pdb

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
    neighs = np.array([[ 1, 0, 0],
                       [-1, 0, 0],
                       [ 0, 1, 0],
                       [ 0,-1, 0],
                       [ 0, 0, 1],
                       [ 0, 0,-1],
                       [ 1, 1, 1],
                       [ 1, 1,-1],
                       [ 1,-1, 1],
                       [ 1,-1,-1],
                       [-1, 1, 1],
                       [-1, 1,-1],
                       [-1,-1, 1],
                       [-1,-1,-1]])

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
        print '-----\nFluid Site Count = ' + str(np.sum(self.Domain.BlockFluidSiteCounts)) + '\n-----\n'
    
        # CONSTRAINT 1: the length of the file must be equal to the
        # value we calculate from the headers.
        if (self.PreambleBytes + self.HeaderBytes + np.sum(self.BlockDataLength)) != os.path.getsize(self.FileName):
            self.PrintError('ERROR: File length appears incorrect.')

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

    def PrintErrors(self, errors):
        for e in errors:
            print 'ERROR: ' + e
            continue
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
        self.errors = []
        self.numFluid = 0
        
        dom = block.Domain
        
        if isinstance(block, AllSolidBlock):
            self.CheckFluidSiteCount()
            return self.errors
        
        i = block.Index
        j = np.zeros(3, dtype=int)

        minInd = i*dom.BlockSize
        maxInd = (i+1)*dom.BlockSize

        for j[0] in xrange(minInd[0], maxInd[0]):
            for j[1] in xrange(minInd[1], maxInd[1]):
                for j[2] in xrange(minInd[2], maxInd[2]):

                    site = block.GetSite(j)
                    type = site.Type
                    edge = site.Edge

                    if not (site.Config == cfg.SOLID_TYPE or site.Config == cfg.FLUID_TYPE):
                        isTypeKnown = False
                        if type == cfg.INLET_TYPE or type == cfg.OUTLET_TYPE:
                            isTypeKnown = True
                            pass

                        if edge:
                            isTypeKnown = True
                            pass

                        if not isTypeKnown:
                            self.errors.append('Site doesn\'t appear to have any fitting type.')
                            pass
                        pass

                    if type == cfg.SOLID_TYPE:
                       continue
                    else:
                        self.numFluid += 1
                        pass

                    nonFluidNeighs = 0
                    offender = np.array([0,0,0])
                    invalid = False
                    for delta in self.neighs:
                        neigh = j + delta
                        if min(neigh) < 0 or min(self.domainExtremeSite - neigh) < 0:
                            nonFluidNeighs += 1
                            offender = neigh
                            invalid = True
                            break

                        elif block.GetSite(neigh).Type == cfg.SOLID_TYPE:
                            nonFluidNeighs += 1
                            offender = neigh
                            break
                        continue

                    # CON 3
                    if type == cfg.FLUID_TYPE and not edge:
                        if nonFluidNeighs != 0:
                            self.errors.append(
                                         'Fluid site at ' + str(j) + ' has type, edge of ' + str(type) + ', ' + str(edge) + ' and has some neighbours that aren\'t valid or are solids, e.g. site at ' + str(offender) + (' which is out-of-bounds' if invalid else ('with type ' + str(block.GetSite(offender).Type))))
                            pass
                        pass

                    # CON 4
                    if type != cfg.FLUID_TYPE or edge:
                        if nonFluidNeighs <= 0:
                            self.errors.append('Site is non-normal fluid but has all valid fluid neighbours.')
                            pass
                        pass

                    if type != cfg.FLUID_TYPE:
                        # CON 5
                        if min(site.CutDistances >= 1.):
                            self.errors.append('Non-fluid site had no cut distances between 0 and 1')

                        for cutDist in site.CutDistances:
                            if cutDist < 0. or (cutDist > 1 and cutDist != float('inf')):
                                # CON 6
                                self.errors.append('Invalid cut distance of ' + str(cutDist))

                    if type == cfg.INLET_TYPE or type == cfg.OUTLET_TYPE:
                        normMag = np.sum(site.BoundaryNormal ** 2) ** 0.5
                        # CON 7
                        if abs(normMag - 1) >= 0.001:
                            self.errors.append('Boundary normal had non-unity magnitude: ' + str(site.BoundaryNormal))
                        # CON 8
                        if site.BoundaryDistance <= 0. or site.BoundaryDistance >= 1.:
                            self.errors.append('Boundary distance was not in [0,1]: ' + str(site.BoundaryDistance))

                    if edge:
                        normMag = np.sum(site.WallNormal ** 2) ** 0.5
                        # CON 9
                        if abs(normMag - 1) >= 0.001:
                            self.errors.append('Wall normal had non-unity magnitude: ' + str(site.WallNormal))
                        # CON 10
                        if site.WallDistance <= 0. or site.WallDistance >= (3.0**0.5):
                            self.errors.append('Wall distance was not in [0,root(3)]: ' + str(site.WallDistance))
                    continue
                continue
            continue
        self.CheckFluidSiteCount()
        return self.errors
    
    def CheckFluidSiteCount(self):
        # CON 2
        if self.numFluid != self.block.nFluidSites:
            self.errors.append('Fluid site counts not equal.')
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


    


