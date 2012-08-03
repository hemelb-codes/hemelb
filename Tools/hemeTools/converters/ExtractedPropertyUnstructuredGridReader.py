"""Define a filter to add the data in a HemeLB extracted property file to the geometry
defined by a HemeLB geometry file.

This module can be run as a script on the command line to convert an extraction to
the corresponding .vtu (VTk Unstructured grid XML file), for usage, run the
script with no arguments. 
"""

import vtk
from hemeTools.parsers.extraction import ExtractedProperty

class ExtractedPropertyUnstructuredGridReader(vtk.vtkProgrammableFilter):
    """VTK-style filter for reading HemeLB extracted property files as VTK data into the
    geometry provided as input. This input must be that derived from the HemeLB
    geometry file used by the simulation from which the extraction file comes and be 
    as output by GmyUnstructuredGridReader (e.g. no scaling), as this class 
    uses the positions to match points.
    
    The vtkUnstructuredGrid will have the same points and cells as the input
    but it will have cell data corresponding to the extraction. The fields take
    the names given in the extracted property file with the same units.
    Position units are metres. The object has no point data.
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
        extraction = ExtractedProperty(self.FileName)
        
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
        grid.ShallowCopy(input)

        extracted_data = extraction.GetByTimeStep(-1)        

        nCells = len(extracted_data.position)

        # Create the field datasets to write into
        field_dict = {}
  
        for read_field in extraction.GetFieldSpec():
            field_name = read_field.name

            field = vtk.vtkDoubleArray()
            field.SetNumberOfComponents(read_field.length)
            field.SetNumberOfTuples(nCells)
            field.SetName(field_name)

            field_dict[field_name]=field

        # The hard work bit.
        for point in extracted_data:
            # Get cell in geometry corresponding to the point
            cellId = locator.FindClosestPoint(point.position)
            if cellId == -1:
                raise ValueError("Can't find cell for point at " + str(point.position))

            for field_name,field in field_dict:
                if field.GetNumberOfComponents() == 3:
                    field.SetTuple3(cellId, *GetAttr(point, field_name))
                else:
                    field.SetTuple1(cellId, GetAttr(point, field_name))
            continue
        
        # Add the arrays to the output
        for field in field_dict:
            grid.GetCellData().AddArray(field_dict[field])
        
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
    parser.add_argument('extracted', nargs=1,
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
    
    converter = ExtractedPropertyUnstructuredGridReader()
    converter.SetInputConnection(reader.GetOutputPort())
    
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputConnection(converter.GetOutputPort())
    
    converter.SetFileName(args.extracted[0])
        
    base, ext = os.path.splitext(args.extracted[0])
    writer.SetFileName(base + '.vtu')
    writer.Write()    
