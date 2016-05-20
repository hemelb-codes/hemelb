import numpy as np
import vtk
import sys
import os.path

def NumpyToVtkImageData(infile):
    array = np.load(infile)

    shape = array.shape

    dataImporter = vtk.vtkImageImport()
    data_matrix = array.astype(np.uint8)
    # Note the transpose to deal with VTK's Fortran style ordering
    data_string = data_matrix.transpose((2,1,0)).tostring()
    dataImporter.CopyImportVoidPointer(data_string, len(data_string))

    dataImporter.SetDataScalarTypeToUnsignedChar()
    dataImporter.SetNumberOfScalarComponents(1)

    extent = (0, shape[0]-1,
              0, shape[1]-1,
              0, shape[2]-1)

    dataImporter.SetDataExtent(*extent)
    dataImporter.SetWholeExtent(*extent)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetInputConnection(dataImporter.GetOutputPort())

    outfile = os.path.splitext(infile)[0] + '.vti'
    writer.SetFileName(outfile)
    writer.Write()
    return

if __name__ == "__main__":
    for infile in sys.argv[1:]:
        NumpyToVtkImageData(infile)
        