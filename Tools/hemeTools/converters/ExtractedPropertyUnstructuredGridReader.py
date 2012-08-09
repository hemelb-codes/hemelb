"""Define a filter to add the data in a HemeLB extracted property file to the geometry
defined by a HemeLB geometry file.

This module can be run as a script on the command line to convert an extraction to
the corresponding .vtu (VTk Unstructured grid XML file), for usage, run the
script with no arguments. 
"""

import vtk
from hemeTools.parsers.extraction import ExtractedProperty

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
        return
    
    def SetExtraction(self, extracted):
        """The parsed extracted property file.
        """
        self.Extracted = extracted
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

    def _Execute(self):
        """Private method that actually does the reading. Called by the VTK
        API.
        """
        input = self.GetUnstructuredGridInput()
        
        # Get the centres as these should match the extracted property positions
        centers = vtk.vtkCellCenters()
        centers.SetInput(input)
        centers.Update()
        # Use this to find the cell ID for each point.
        locator = vtk.vtkOctreePointLocator()
        locator.SetDataSet(centers.GetOutput())
        locator.BuildLocator()
        # Should be fine enough
        locator.SetTolerance(0.1 * extraction.voxelSizeMetres)
        
        # Copy the structure to output.
        grid = self.GetUnstructuredGridOutput()

        # Get the data we want to visualise, Note that this is likely to 
        # involve a subset of the points in the geometry
        extracted_data = self.Extracted.GetByTimeStep(self.Time)

        # Make a list of the cell ids to keep
        cellIds = vtk.vtkIdTypeArray()
        cellIds.SetNumberOfComponents(1)
 
        for point in extracted_data:
            # Get cell in geometry corresponding to the point
            cellId = locator.FindClosestPoint(point.position)

            if cellId == -1:
                raise ValueError("Can't find cell for point at " + str(point.position))

            cellIds.InsertNextValue(cellId)

        # Make an object to select only the cell ids we want
        selector = vtk.vtkSelectionNode()
        selector.SetFieldType(vtk.vtkSelectionNode.CELL)
        selector.SetContentType(vtk.vtkSelectionNode.INDICES)
        selector.SetSelectionList(cellIds)

        # Make an object to hold the selector
        selectors = vtk.vtkSelection()
        selectors.AddNode(selector)

        # Perform the selection
        extractSelection = vtk.vtkExtractSelectedIds()
        extractSelection.SetInput(0, input)
        extractSelection.SetInput(1, selectors)
        extractSelection.Update()

        # Copy the result into our grid, and get the lookup array from the new cell ids
        # to the old cell ids.
        grid.ShallowCopy(extractSelection.GetOutput())

        originalCellIdLookup = grid.GetCellData().GetArray('vtkOriginalCellIds')

        nCells = grid.GetNumberOfCells()

        # Create the field datasets to write into
        field_dict = {}
  
        for name, xdrType, memType, length, offset in self.Extracted.GetFieldSpec():
            # Skip the grid.
            if name == 'grid':
                continue

            # Create a VTK object for storing the data.
            field = vtk.vtkDoubleArray()
            if isinstance(length, tuple):
                length = length[0]
            field.SetNumberOfComponents(length)
            field.SetNumberOfTuples(nCells)
            field.SetName(name)

            # Insert it into the dictionary
            field_dict[name] = field

        # The hard work bit.
        for point in extracted_data:
            # Get cell in geometry corresponding to the point
            oldCellId = locator.FindClosestPoint(point.position)
            if oldCellId == -1:
                raise ValueError("Can't find cell for point at " + str(point.position))

            cellId = -1

            # Translate to new grid.
            for i in range(0, originalCellIdLookup.GetSize()):
                if originalCellIdLookup.GetComponent(0, i) == oldCellId:
                    cellId = i
                    break

            if cellId == -1:
                raise ValueError("Can't translate old cellId %i to new grid for point at %s" % (oldCellId, str(point.position)))

            # For each field considered, set the property in the associated VTK array.
            for field_name,field in field_dict.items():
                if field.GetNumberOfComponents() == 3:
                    field.SetTuple3(cellId, *getattr(point, field_name))
                else:
                    field.SetTuple1(cellId, getattr(point, field_name))
        
        # Add the arrays to the output
        for field in field_dict:
            grid.GetCellData().AddArray(field_dict[field])

        # Remove the arrays added during the translation of the grid.
        grid.GetCellData().RemoveArray('vtkOriginalCellIds')
        grid.GetPointData().RemoveArray('vtkOriginalPointIds')

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
        extraction=ExtractedProperty(extractedFile)

        # For each timestamp in that file, create a reader...
        for time in extraction.times:
            converter = ExtractedPropertyUnstructuredGridReader()
            converter.SetInputConnection(reader.GetOutputPort())
            converter.SetExtraction(extraction)
            converter.SetTime(time)
    
            # And a writer...
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetInputConnection(converter.GetOutputPort())
        
            # And give the output a unique filename based on the extracted
            # property file name and the timestamp.
            base, ext = os.path.splitext(extractedFile)
            writer.SetFileName('%s_%s.vtu' % (base, str(time)))
            writer.Write()
