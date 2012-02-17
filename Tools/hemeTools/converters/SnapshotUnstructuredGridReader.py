import vtk
from hemeTools.parsers.snapshot import HemeLbSnapshot

class SnapshotUnstructuredGridReader(vtk.vtkProgrammableFilter):
    def __init__(self):
        self.FileName = ""
        self.SetExecuteMethod(self.Execute)
        
        return
    def SetFileName(self, path):
        self.FileName = path
        return
    def GetFileName(self):
        return self.FileName
    
    def Execute(self):
        input = self.GetUnstructuredGridInput()
         
        centers = vtk.vtkCellCenters()
        centers.SetInput(input)
        centers.Update()
        
        locator = vtk.vtkOctreePointLocator()
        locator.SetDataSet(centers.GetOutput())
        locator.BuildLocator()
        
        locator.SetTolerance(0.1 * snap.voxel_size)
        
        snap = HemeLbSnapshot(self.FileName)
        
        grid = self.GetUnstructuredGridOutput()
        grid.ShallowCopy(input)
        
        nCells = len(snap)
        if grid.GetNumberOfCells() != nCells:
            raise ValueError('Skeleton ({}) and snapshot ({}) have different cell counts'.format(grid.GetNumberOfCells(),
                                                                                                 len(snap)))
        velocity = vtk.vtkDoubleArray()
        velocity.SetNumberOfComponents(3)
        velocity.SetNumberOfTuples(nCells)
        velocity.SetName('velocity')
        
        pressure = vtk.vtkDoubleArray()
        pressure.SetNumberOfComponents(1)
        pressure.SetNumberOfTuples(nCells)
        pressure.SetName('pressure')
        
        stress = vtk.vtkDoubleArray()
        stress.SetNumberOfComponents(1)
        stress.SetNumberOfTuples(nCells)
        stress.SetName('stress')
        
        for pt in snap:
            cellId = locator.FindClosestPoint(pt.position)
            if cellId == -1:
                raise ValueError("Can't find cell for point at " + str(pt.position))
            
            velocity.SetTuple3(cellId, *pt.velocity)
            pressure.SetTuple1(cellId, pt.pressure)
            stress.SetTuple1(cellId, pt.stress)
            continue
        
        grid.GetCellData().AddArray(velocity)
        grid.GetCellData().AddArray(pressure)
        grid.GetCellData().AddArray(stress)
        return grid
    pass

if __name__ == '__main__':
    import sys
    import os.path
    """At the moment this works like:
    ./script.py skelton.vtu snap1.dat snap2.dat ...
    
    This should remain.
    However the code below must be made similar to GmyUGReader.py
    Must connect a vtkAlgorithm to the input of above filter
    either GmyUGReader or a vtkXMLUnstrGridReader
    and then set the filename in the for loop.
    """
    
    f = Flesher(sys.argv[1])
    writer = vtk.vtkXMLUnstructuredGridWriter()
    for snapName in sys.argv[2:]:
        base, ext = os.path.splitext(snapName)
        writer.SetFileName(base + '.vtu')

        snap = HemeLbSnapshot(snapName)
        ug = f(snap)
        writer.SetInput(ug)
        writer.Write()

        
