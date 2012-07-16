# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

"""This module contains a vtkProgrammableSource subclass that will read a
HemeLB Geometry file (.gmy) and create a VTK data structure with the same
geometry. 

This module can be run as a script on the command line to convert a .gmy file
to the corresponding .vtu (VTk Unstructured grid XML file), for usage, run the
script with no arguments.

"""
import vtk
import numpy as np
from hemeTools.parsers.geometry.simple import ConfigLoader

class GmyUnstructuredGridReader(vtk.vtkProgrammableSource):
    """VTK-style reader for HemeLB Geometry files (.gmy). When run, it will 
    create a VTK data structure with the same geometry as the file.
    
    The vtkUnstructuredGrid that is the output will have one cell for every 
    fluid site in the geometry file. The cell will be a cube or voxel (in VTK 
    terminology a cell with type VTK_VOXEL or, when accessed through the 
    GetCell() API, an instance of vtkVoxel) with its centre on the lattice 
    point and its corners half a lattice unit away in each of the eight 
    <1,1,1> type directions.
    
    The output units are metres and the object has no cell data or point data.
    """
    def __init__(self):
        self.SetExecuteMethod(self._Execute)
        self.GetUnstructuredGridOutput()
        self.FileName = ""
        return
    
    def GetOutputPort(self, index=3):
        return vtk.vtkAlgorithm.GetOutputPort(self, index)
    
    def SetFileName(self, name):
        """The file to read.
        """
        self.FileName = name
        self.Modified()
        return
    
    def GetFileName(self):
        """The file to read.
        """
        return self.FileName

    def _Execute(self):
        """Private method that actually does the reading. Called by the VTK
        API.
        """
        output = self.GetUnstructuredGridOutput()
        loader = GeneratingLoader(self.FileName)
        loader.Load()
        
        output.ShallowCopy(loader.Grid)
        return
    
class GeneratingLoader(ConfigLoader):
    """This subclass will create a vtkUnstructuredGrid of all the fluid sites
    as is loads in a geometry file.
    
    It is intended only for use by the GmyUnstructedGridReader. However, for 
    reference, after it has loaded the geometry file (through a call to Load(),
    the vtkUnstructuredGrid will be in the "Grid" attribute.
    
    """
        
    def OnEndHeader(self):
        """Get ready to add blocks
        """
        self.TotalFluidSites = self.Domain.BlockFluidSiteCounts.sum()
        
        # Efficiently add the voxel corner points
        self.Locator = vtk.vtkMergePoints()
        self.Points = vtk.vtkPoints()
        
        # Compute the bounds of the voxel corner points
        # Initially in a shape == (3,2) array [[0, nSites[x]], ... ]
        bounds = np.concatenate((np.zeros(3,dtype=int),
                                 self.Domain.SiteCounts)).reshape((2,3)).transpose()
        # shift half a voxel down, scale and then translate
        bounds = self.Domain.VoxelSize*(bounds - 0.5) + \
            self.Domain.Origin[:, np.newaxis]
        
        # Since we know how the points will be arranged, use that information
        # to tune the vtkMergePoints
        self.Locator.SetDivisions(*self.Domain.BlockCounts)
        self.Locator.InitPointInsertion(self.Points, bounds.flatten(),
                                        self.TotalFluidSites)
        
        # This is what we're trying to construct. Allocate data (we know how
        # many cells (= TotalFluidSites) but not how many points, as we don't
        # know the number of surface sites.
        self.Grid = vtk.vtkUnstructuredGrid()
        self.Grid.Allocate(self.TotalFluidSites, 1000)
        self.Grid.SetPoints(self.Points)
        
        # Array of offsets from a point to the corners of its cell.
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
        # Used as a working array.
        self.vox = vtk.vtkVoxel()
        
        return

    def OnEndSite(self, block, site):
        """Called once a site has been read from the file.
        """
        if site.IsSolid:
            return
        
        # Save two attribute accesses and a function call
        idlist = self.vox.GetPointIds()
        
        for i, delta in enumerate(self.deltas):
            # For each of the corner points, see if it's already added
            # If so, just use its ID, if not insert it into the array and 
            # get the new ID.
            pt = site.Position + delta
            ptId = self.Locator.IsInsertedPoint(pt)
            if ptId == -1:
                ptId = self.Locator.InsertNextPoint(pt)
                pass
            
            idlist.SetId(i, ptId)
            continue
        
        # Insert the cell into the unstructured grid.
        self.Grid.InsertNextCell(self.vox.GetCellType(), idlist)
        return
    
    def OnEndBlock(self, bIdx, bIjk):
        """Since we're done with the block, we can free its sites.
        """
        self.Domain.GetBlock(bIdx).DeleteSites()
        return
    
    def OnEndBody(self):
        """We are done, so free all we can.
        """
        grid = self.Grid
        self.__dict__.clear()
        self.Grid = grid
        return

    
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Create a VTK Unstructured grid with the correct topology '
        'for the given geometry file.'
        )
    parser.add_argument('input', nargs=1, type=argparse.FileType(),
                        help='the geometry file')
    parser.add_argument('output', nargs=argparse.OPTIONAL,
                        help='file to write VTK data to, if omitted, use '
                        'input basename + ".vtu"',
                        default=None)
    
    args = parser.parse_args()
    
    import os.path
    # We have to have the file name, not the open file, but argparse will give
    # nicer usage instructions if we tell it to open for us. Replace the file
    # handle with its name.
    input = args.input[0].name
    args.input[0].close()
    
    # If output wasn't given, use the input with '.vtu' extension.
    if args.output is None:
        base, ext = os.path.splitext(input)
        output = base + '.vtu'
    else:
        output = args.output
        pass
    
    reader = GmyUnstructuredGridReader()
    reader.SetFileName(input)
    
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(output)
    
    writer.SetInputConnection(reader.GetOutputPort())
    
    writer.Write()
    
