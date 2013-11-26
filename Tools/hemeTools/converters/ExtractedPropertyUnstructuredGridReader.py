"""Define a filter to add the data in a HemeLB extracted property file to the geometry
defined by a HemeLB geometry file.

This module can be run as a script on the command line to convert an extraction to
the corresponding .vtu (VTk Unstructured grid XML file), for usage, run the
script with no arguments. 
"""

import vtk
from vtk.util import numpy_support
import numpy as np
from hemeTools.utils import MatchCorresponding
from hemeTools.parsers.extraction import ExtractedProperty

# Work around limitations of the VTK bindings support (enum wrapping added in 5.8)
v = vtk.vtkVersion()
if v.GetVTKMajorVersion() >5 or \
    (v.GetVTKMajorVersion()==5 and v.GetVTKMinorVersion() >6):
    CELL = vtk.vtkSelectionNode.CELL
    INDICES = vtk.vtkSelectionNode.INDICES
else:
    CELL = 0
    INDICES = 4
del v

class ExtractedPropertyUnstructuredGridReader(vtk.vtkProgrammableFilter):
    """VTK-style filter for reading HemeLB extracted property files as VTK 
    data into the geometry provided as input. This input must be that derived 
    from the HemeLB geometry file used by the simulation from which the 
    extraction file comes and be as output by GmyUnstructuredGridReader (e.g. 
    no scaling), as this class uses the positions to match points.
    
    The vtkUnstructuredGrid will have that subset of the points and cells from
    the input that contain data. Cell data will be as in the extraction file. 
    The fields take the names and units given in the extracted property file.
    Position units are metres. The object has no point data.
    """
    def __init__(self):
        self.SetExecuteMethod(self._Execute)
        self.Extracted = None
        self.Time = None
        self.Skeleton = None
        self.OutputCellIdsByInputIndex = None
        self.InputMTime = -1
        return
    
    def SetExtraction(self, extracted):
        """The parsed extracted property file.
        """
        self.Extracted = extracted
        self.Skeleton = None
        self.OutputCellIdsByInputIndex = None
        self.Modified() 
        return

    def SetTime(self, time):
        """The timestamp to read properties for.
        """
        time = int(time)
        if self.Time != time:
            self.Time = time
            self.Modified()
        return    

    def _CreateSkeleton(self):
        """Create the structure of the output vtkUnstructuredGrid and a map
        from the index of a point in the extraction file to the corresponding
        cellId in the skeleton.
        
        This method should only be called if the extraction object this 
        instance is working on has changed since the last time this method was
        called.
        """
        input = self.GetUnstructuredGridInput()
        
        # Get the centres as these should match the extracted property positions
        centers = vtk.vtkCellCenters()
        if vtk.vtkVersion().GetVTKMajorVersion() <= 5:
          centers.SetInput( input );
        else:
          centers.SetInputData( input );
        #centers.SetInput(input)
        centers.Update()
        # Use this to find the cell ID for each point.
        locator = vtk.vtkOctreePointLocator()
        locator.SetDataSet(centers.GetOutput())
        locator.BuildLocator()
        # Should be fine enough
        locator.SetTolerance(0.1 * self.Extracted.voxelSizeMetres)

        # Get the first set of extracted data from the file. We don't care
        # which as we only want to use the positions.
        extracted_positions = self.Extracted.GetByIndex(0).position
        nExtractedPoints = len(extracted_positions)
        
        # Make a list of the cell ids to keep; i.e. the cell in the input 
        # (whole geometry) with ID cellIdsGmy[i] contains  the point given by
        # extracted_positions[i]
        gmyCellIdsByInput = vtk.vtkIdTypeArray()
        gmyCellIdsByInput.SetNumberOfComponents(1)
        gmyCellIdsByInput.SetNumberOfTuples(nExtractedPoints)
        
        
        gmyCellIdsByInputNp = numpy_support.vtk_to_numpy(gmyCellIdsByInput)
         
        for i, point in enumerate(extracted_positions):
            # Get cell in geometry corresponding to the point
            cellId = locator.FindClosestPoint(point)

            if cellId == -1:
                raise ValueError("Can't find cell for point at " + str(point))

            gmyCellIdsByInputNp[i] = cellId
        
        # Make an object to select only the cell ids we want
        selector = vtk.vtkSelectionNode()
        selector.SetFieldType(CELL)
        selector.SetContentType(INDICES)
        selector.SetSelectionList(gmyCellIdsByInput)
        
        # Make an object to hold the selector
        selectors = vtk.vtkSelection()
        selectors.AddNode(selector)
        
        # Perform the selection
        extractSelection = vtk.vtkExtractSelectedIds()
        if vtk.vtkVersion().GetVTKMajorVersion() <= 5:        
            extractSelection.SetInput(0, input)
            extractSelection.SetInput(1, selectors)
        else:
            extractSelection.SetInputData(0, input)
            extractSelection.SetInputData(1, selectors)
        extractSelection.Update()
        
        gmyCellIdsByOutput = extractSelection.GetOutput().GetCellData().GetArray('vtkOriginalCellIds')
        
        
        gmyCellIdsByOutputNp = numpy_support.vtk_to_numpy(gmyCellIdsByOutput)
        
        self.OutputCellIdsByInputIndex = MatchCorresponding(gmyCellIdsByOutputNp,gmyCellIdsByInputNp)
        
        # Create the skeleton data object and only copy in the structure.
        self.Skeleton = vtk.vtkUnstructuredGrid()
        self.Skeleton.CopyStructure(extractSelection.GetOutput())
        return
    
    def _Execute(self):
        """Private method that actually does the reading. Called by the VTK
        API.
        """
        # Check if we have a skeleton and if so, if it is up-to-date.
        if self.Skeleton is None or (self.GetInput().GetMTime() > self.InputMTime):
            self._CreateSkeleton()
            self.InputMTime = self.GetInput().GetMTime()
        
        # Copy the structure to output.
        output = self.GetUnstructuredGridOutput()
        output.CopyStructure(self.Skeleton)
        
        # Get the data we want to visualise, Note that this is likely to 
        # involve a subset of the points in the geometry
        extracted_data = self.Extracted.GetByTimeStep(self.Time)
        nExtractedPoints = len(extracted_data)
        
        # Create the arrays to store field data
        field_dict = {}
        
        for name, xdrType, memType, length, offset in self.Extracted.GetFieldSpec():
            # Skip the grid (i.e. 3d-index)
            if name == 'grid':
                continue
            
            # Set up the VTK array
            field = vtk.vtkDoubleArray()
            if isinstance(length, tuple):
                length = length[0]
            field.SetNumberOfComponents(length)
            field.SetNumberOfTuples(nExtractedPoints)
            field.SetName(name)

            # Insert it into the dictionary
            field_dict[name] = field

        # Copy the data into the correct position in the output array and add 
        # the array to the output.
        # TODO: this needs a case for the stress. We should check what 
        # representation VTK uses for rank 3 tensors. 
        for field_name, field in field_dict.iteritems():
            fieldArray = numpy_support.vtk_to_numpy(field)
            fieldArray[self.OutputCellIdsByInputIndex] = getattr(extracted_data, field_name)
            output.GetCellData().AddArray(field)
                        
        return
    pass

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Create a VTK Unstructured grid with geometry of the '
        'first argument and the data in the extracted property file in all later arguments.'
        )
    parser.add_argument('geometry', nargs=1,
                        help='the geometry, either a HemeLB .gmy file '
                        'or a derived VTK unstructured grid')
    parser.add_argument('extracted', nargs=argparse.ONE_OR_MORE,
                        help='extracted property file to convert to VTK, output will'
                        ' be put in the input basename + ".vtu"')
    
    args = parser.parse_args()
    geometry = args.geometry[0]
    import os.path
    
    base, ext = os.path.splitext(geometry)
    if ext == '.gmy':
        from hemeTools.converters.GmyUnstructuredGridReader import GmyUnstructuredGridReader
        reader = GmyUnstructuredGridReader()
    elif ext == '.vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
         print ('Cannot infer reader time from file extension "%s"' % ext)
         parser.print_help()
         raise SystemExit()
    
    reader.SetFileName(geometry)

    # For each extracted file given, parse the file...    
    for extractedFile in args.extracted:
        extraction = ExtractedProperty(extractedFile)
        converter = ExtractedPropertyUnstructuredGridReader()
        converter.SetInputConnection(reader.GetOutputPort())
        converter.SetExtraction(extraction)
            
        # Create a writer...
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputConnection(converter.GetOutputPort())
        
        # For each timestamp in that file
        for time in extraction.times:
            # Choose the time
            converter.SetTime(time)
            
            # And give the output a unique filename based on the extracted
            # property file name and the timestamp.
            base, ext = os.path.splitext(extractedFile)
            writer.SetFileName('%s_%s.vtu' % (base, str(time)))
            # trigger the update
            writer.Write()
