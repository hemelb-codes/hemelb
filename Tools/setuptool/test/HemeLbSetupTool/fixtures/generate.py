import vtk
from .load import fixture_path

def write_stl(input,label):
    writer = vtk.vtkSTLWriter()
    writer.SetInputConnection(input.GetOutputPort())
    writer.SetFileName(fixture_path(label))
    writer.Write()

def cube():
    source=vtk.vtkCubeSource()
    triangulator=vtk.vtkTriangleFilter()
    triangulator.SetInputConnection(source.GetOutputPort())
    # Make some polydata
    write_stl(triangulator,"cube")

