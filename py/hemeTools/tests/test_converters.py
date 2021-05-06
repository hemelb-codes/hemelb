import os.path

import numpy as np
from vtk import vtkUnstructuredGrid, vtkVoxel

# from hemeTools.converters import ExtractedPropertyUnstructuredGridReader as xtr
from hemeTools.converters import GmyUnstructuredGridReader as gmy


def check_difftest_gmy(grid):
    assert isinstance(grid, vtkUnstructuredGrid)
    # No data
    assert grid.GetCellData().GetNumberOfArrays() == 0
    assert grid.GetFieldData().GetNumberOfArrays() == 0
    assert grid.GetPointData().GetNumberOfArrays() == 0

    n_sites = 44250
    assert grid.GetNumberOfCells() == n_sites
    for i in range(n_sites):
        s = grid.GetCell(i)
        assert isinstance(s, vtkVoxel)
        # Would be nice to check inside the cylinder but skip for now


def test_converting_gmy_loader(diffTestDir):
    config_path = os.path.join(diffTestDir, "config.xml")
    loader = gmy.GeneratingLoader(config_path)
    loader.Load()
    check_difftest_gmy(loader.Grid)


def test_GmyUnstructuredGridReader(diffTestDir):
    config_path = os.path.join(diffTestDir, "config.xml")
    reader = gmy.GmyUnstructuredGridReader()
    reader.SetFileName(config_path)
    reader.Update()
    check_difftest_gmy(reader.GetOutputDataObject(0))
