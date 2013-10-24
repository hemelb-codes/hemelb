from vtk.util.numpy_support import vtk_to_numpy

def IterCellPtIds(cellArray):
    """Iterate over the cells, yielding NumPy arrays of the point IDs
    making up each cell.
    """
    data = vtk_to_numpy(cellArray.GetData())
    n = len(data)
    
    lenId = 0
    while lenId < n:
        start = lenId + 1
        end = start + data[lenId]
        
        yield data[start:end]

        lenId = end
        continue
    return
