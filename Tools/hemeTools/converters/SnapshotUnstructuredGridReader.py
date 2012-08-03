# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

"""Define a filter to add the data in a HemeLB snapshot file to the geometry
defined by a HemeLB geometry file.

This module can be run as a script on the command line to convert snapshot to
the corresponding .vtu (VTk Unstructured grid XML file), for usage, run the
script with no arguments. 
"""

import vtk
from hemeTools.parsers.snapshot import HemeLbSnapshot

class SnapshotUnstructuredGridReader(vtk.vtkProgrammableFilter):
    """VTK-style filter for reading HemeLB snapshot files as VTK data into the
    geometry provided as input. This input must be that derived from the HemeLB
    geometry file used by the simulation from which the snapshot comes an be 
    as output by GmyUnstructuredGridReader (e.g. no scaling), as this class 
    uses the positions to match points.
    
    The vtkUnstructuredGrid will have the same points and cells as the input
    but it will have cell data corresponding to the snapshot. The fields are
    named 'velocity', 'pressure' and 'stress' with units of (respectively)
    m/s, mmHg and Pa. Position units are metres. The object has no point data.
    """
    def __init__(self):
        self.FileName = ""
        self.SetExecuteMethod(self._Execute)
        
        return
    
    def SetFileName(self, path):
        """The file to read.
        """
        self.FileName = path
        return
    
    def GetFileName(self):
        """The file to read.
        """
        return self.FileName
    
    def _Execute(self):
        """Private method that actually does the reading. Called by the VTK
        API.
        """
        input = self.GetUnstructuredGridInput()
        
        # Load the snapshot data
        snap = HemeLbSnapshot(self.FileName)
        
        # Get the centres as these should match the snapshot positions
        centers = vtk.vtkCellCenters()
        centers.SetInput(input)
        centers.Update()
        # Use this to find the cell ID for each point.
        locator = vtk.vtkOctreePointLocator()
        locator.SetDataSet(centers.GetOutput())
        locator.BuildLocator()
        # Should be fine enough
        locator.SetTolerance(0.1 * snap.voxel_size)
        
        # Copy the structure to output.
        grid = self.GetUnstructuredGridOutput()
        grid.ShallowCopy(input)
        
        nCells = len(snap)
        # Basic sanity check that we probably have matching geometry and
        # snapshot.
        if grid.GetNumberOfCells() != nCells:
            raise ValueError('Geometry ({}) and snapshot ({}) have different '
                             'cell counts'.format(grid.GetNumberOfCells(),
                                                  len(snap)))
        # Velocity field
        velocity = vtk.vtkDoubleArray()
        velocity.SetNumberOfComponents(3)
        velocity.SetNumberOfTuples(nCells)
        velocity.SetName('velocity')
        
        # Pressure
        pressure = vtk.vtkDoubleArray()
        pressure.SetNumberOfComponents(1)
        pressure.SetNumberOfTuples(nCells)
        pressure.SetName('pressure')
        
        # Stress
        stress = vtk.vtkDoubleArray()
        stress.SetNumberOfComponents(1)
        stress.SetNumberOfTuples(nCells)
        stress.SetName('stress')
        
        # The hard work bit.
        for pt in snap:
            # Get cell in geometry corresponding to the point
            cellId = locator.FindClosestPoint(pt.position)
            if cellId == -1:
                raise ValueError("Can't find cell for point at " + str(pt.position))
            
            # Copy the data into the VTK structure
            velocity.SetTuple3(cellId, *pt.velocity)
            pressure.SetTuple1(cellId, pt.pressure)
            stress.SetTuple1(cellId, pt.stress)
            continue
        
        # Add the arrays to the output
        grid.GetCellData().AddArray(velocity)
        grid.GetCellData().AddArray(pressure)
        grid.GetCellData().AddArray(stress)
        
        return
    pass

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Create a VTK Unstructured grid with geometry of the '
        'first argument and the data in the snapshots in all later arguments.'
        )
    parser.add_argument('geometry', nargs=1,
                        help='the geometry, either a HemeLB .gmy file '
                        'or a derived VTK unstructured grid')
    parser.add_argument('snapshots', nargs=argparse.ONE_OR_MORE,
                        help='snapshot file(s) to convert to VTK, output will'
                        ' be put in the input basename + ".vtu"')
    
    args = parser.parse_args()
    geometry = args.geometry[0]
    import os.path
    
    base, ext = os.path.splitext(geometry)
    if ext == '.gmy':
        from .GmyUnstructuredGridReader import GmyUnstructuredGridReader
        reader = GmyUnstructuredGridReader()
    elif ext == '.vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
         print ('Cannot infer reader time from file extension "%s"' % ext)
         parser.print_help()
         raise SystemExit()
    
    reader.SetFileName(geometry)
    
    converter = SnapshotUnstructuredGridReader()
    converter.SetInputConnection(reader.GetOutputPort())
    
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputConnection(converter.GetOutputPort())
    
    for snapName in args.snapshots:
        converter.SetFileName(snapName)
        
        base, ext = os.path.splitext(snapName)
        writer.SetFileName(base + '.vtu')
        writer.Write()
        continue
    