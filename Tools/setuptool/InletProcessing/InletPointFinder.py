import numpy as np
import vtk
from hemeTools.parsers.geometry.simple import ConfigLoader
from hemeTools.parsers.geometry.generic import Site

class GeometryInletPointFinder(ConfigLoader):
    """Class to find all the INLET points contained within a geometry file 
    (.gmy). Data is returned as a vtkPolyData, with the 3D index of the points
    supplied as CellData in an array named "PointsIds".
    
    It's based on the hemeTools.parsers.geometry classes.
    """
    
    def __init__(self, filename):
        """filename = name of geometry file to examine
        """
        ConfigLoader.__init__(self, filename)
        self.InletData = {}
        self.InletPointIndices = {}
        self.InletPointPositions = {}
        
        self.Load()
        
    def GetInletData(self):
        ids = self.InletPointIndices.keys()
        ids.sort()
        ans = []
        for i, inletId in enumerate(ids):
            assert i == inletId
            pd = vtk.vtkPolyData()
            pd.SetPoints(self.InletPointPositions[inletId])
            
            self.InletPointIndices[inletId].SetName("PointIds")
            pd.GetPointData().AddArray(self.InletPointIndices[inletId])
            ans.append(pd)
            continue
        
        return ans
    
    @staticmethod
    def _MakeArray(input, arrayType):
        output = arrayType()
        output.SetNumberOfComponents(3)
        output.SetNumberOfTuples(len(input))
        for i, line in enumerate(input):
            output.SetTuple3(i, line)
        return output
    
    def _AddPointForInletId(self, inletId, site):
        """Private method. Adds the point specified by 'site' to the output 
        for the inlet with ID 'inletId'
        """
        try:
            ids = self.InletPointIndices[inletId]
            pos = self.InletPointPositions[inletId]
        except KeyError:
            ids = self.InletPointIndices[inletId] = vtk.vtkIntArray()
            ids.SetNumberOfComponents(3)
            pos = self.InletPointPositions[inletId] = vtk.vtkPoints()
            pass
        
        #ids.append(site.Index)
        #pos.append(site.Position)
        ids.InsertNextTuple3(*site.Index)
        pos.InsertNextPoint(*site.Position)
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
    