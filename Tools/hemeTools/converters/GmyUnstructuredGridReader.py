from vtk import (vtkMergePoints, vtkPoints, vtkUnstructuredGrid, vtkVoxel, 
                 vtkProgrammableSource, vtkXMLUnstructuredGridWriter) 
import numpy as np
from hemeTools.parsers.geometry.simple import ConfigLoader

class GeneratingLoader(ConfigLoader):
    def __init__(self, filename):
        ConfigLoader.__init__(self, filename)
        return
    
    def OnEndHeader(self):
        self.TotalFluidSites = self.Domain.BlockFluidSiteCounts.sum()

        self.Locator = vtkMergePoints()
        self.Points = vtkPoints()
        bounds = np.concatenate((np.zeros(3,dtype=int),
                                 self.Domain.SiteCounts)).reshape((2,3)).transpose()
        bounds = self.Domain.VoxelSize*(bounds - 0.5) + self.Domain.Origin[:, np.newaxis]
        
        self.Locator.SetDivisions(*self.Domain.BlockCounts)
        self.Locator.InitPointInsertion(self.Points, bounds.flatten(), self.TotalFluidSites)
        
        self.Grid = vtkUnstructuredGrid()
        self.Grid.Allocate(self.TotalFluidSites, 1000)
        self.Grid.SetPoints(self.Points)

        self.deltas = 0.5 * self.Domain.VoxelSize * \
            np.array([[ 1, 1, 1],
                      [ 1, 1,-1],
                      [ 1,-1, 1],
                      [ 1,-1,-1],
                      [-1, 1, 1],
                      [-1, 1,-1],
                      [-1,-1, 1],
                      [-1,-1,-1]],
                dtype=float)
        self.vox = vtkVoxel()
        
        return

    def OnEndSite(self, block, site):
        if site.IsSolid:
            return
        
        idlist = self.vox.GetPointIds()
        for i, delta in enumerate(self.deltas):
            pt = site.Position + delta
            ptId = self.Locator.IsInsertedPoint(pt)
            if ptId == -1:
                ptId = self.Locator.InsertNextPoint(pt)
                pass
            idlist.SetId(i, ptId)
            continue
        
        self.Grid.InsertNextCell(self.vox.GetCellType(), idlist)
        return
    
    def OnEndBlock(self, bIdx, bIjk):
        self.Domain.GetBlock(bIdx).DeleteSites()
        return
    
    def OnEndBody(self):
        del self.Domain
        return

class GmyUnstructuredGridReader(vtkProgrammableSource):
    def __init__(self):
        self.SetExecuteMethod(self.Execute)
        self.FileName = ""
        return
    
    def SetFileName(self, name):
        self.FileName = name
        return
    def GetFileName(self):
        return self.FileName

    def Execute(self):
        output = self.GetUnstructuredGridOutput()
        loader = GeneratingLoader(self.FileName)
        loader.Load()
        output.ShallowCopy(loader.Grid)
        return
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Create a VTK Unstructured grid with the correct topology for the given geometry file.'
        )
    parser.add_argument('input', nargs=1,
                        help='the geometry file')
    parser.add_argument('output', nargs=1,
                        help='file to write VTK data to')
    
    args = parser.parse_args()

    reader = GmyUnstructuredGridReader()
    reader.SetFileName(input)
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetFileName(args.output)
    writer.SetInputConnection(reader.GetOutputPort())
    
    writer.Write()
    
