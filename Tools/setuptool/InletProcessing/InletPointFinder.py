# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import numpy as np

from hemeTools.parsers.geometry.simple import ConfigLoader
from hemeTools.parsers.geometry.generic import Site

class InletPointFinder(ConfigLoader):
    """Class to find all the INLET points contained within a geometry file 
    (.gmy).
    
    It's based on the hemeTools.parsers.geometry classes.
    """
    
    def __init__(self, filename):
        """filename = name of geometry file to examine
        """
        ConfigLoader.__init__(self, filename)
        self.InletPointIndices = {}
        self.InletPointPositions = {}
        
    def Run(self):
        self.Load()
        ids = self.InletPointIndices.keys()
        ids.sort()
        ans = []
        for i, inletId in enumerate(ids):
            assert i == inletId
            ids = np.array(self.InletPointIndices[inletId])
            pos = np.array(self.InletPointPositions[inletId])
            ans.append((ids, pos))
            continue
        
        return ans
        
    def _AddPointForInletId(self, inletId, site):
        """Private method. Adds the point specified by 'site' to the output 
        for the inlet with ID 'inletId'
        """
        try:
            ids = self.InletPointIndices[inletId]
            pos = self.InletPointPositions[inletId]
        except KeyError:
            ids = self.InletPointIndices[inletId] = []
            pos = self.InletPointPositions[inletId] = []
            pass
        
        ids.append(site.Index)
        pos.append(site.Position)
        return
    
    def OnEndSite(self, block, site):
        """Triggered once a site has been loaded. Check if it's an inlet and, 
        if so, add data to the output.
        """
        if site.IsSolid or not np.any(site.IntersectionType == Site.INLET_INTERSECTION):
            return
        
        inletIds = site.IOletIndex[np.where(site.IntersectionType == Site.INLET_INTERSECTION)]
        
        # We assume here that a lattice site is adjacent to at most one inlet
        assert np.all(inletIds == inletIds[0])
        inletId = inletIds[0]
        
        self._AddPointForInletId(inletId, site)
        return
    
    def OnEndBlock(self, bIdx, bIjk):
        """Triggered at the end of parsing a block. We can delete the block as
        it's no longer needed.
        """
        self.Domain.DeleteBlock(bIdx)
        return
    