from vmtk.vmtksurfaceviewer import vmtkSurfaceViewer
from vmtk.vmtksurfacereader import vmtkSurfaceReader
import vtk
import os.path

def view(surf):
    if isinstance(surf, vtk.vtkPolyData):
        view_vtkPolyData(surf)
    elif isinstance(surf, str):
        if os.path.exists(surf):
            view_file(surf)
            pass
    else:
        raise TypeError("Dinnae ken hae ta voo't!")
    return

def view_vtkPolyData(surf):
    viewer = vmtkSurfaceViewer()
    viewer.Surface = surf
    viewer.Execute()
    return

def view_file(fn):
    reader = vmtkSurfaceReader()
    reader.InputFileName = fn
    reader.Execute()
    return view_vtkPolyData(reader.Output)

    
