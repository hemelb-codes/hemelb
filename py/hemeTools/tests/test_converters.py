import os.path

import numpy as np
from vtk import vtkUnstructuredGrid, vtkVoxel

from hemeTools.converters import ExtractedPropertyUnstructuredGridReader as xtr
from hemeTools.converters import GmyUnstructuredGridReader as gmy
from hemeTools.parsers.extraction import ExtractedProperty


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


def test_ExtractedPropertyUnstructuredGridReader(diffTestDir):
    config_path = os.path.join(diffTestDir, "config.xml")
    gmy_reader = gmy.GmyUnstructuredGridReader()
    gmy_reader.SetFileName(config_path)

    xtr_path = os.path.join(diffTestDir, "CleanExtracted", "flow_snapshot.dat")
    extraction = ExtractedProperty(xtr_path)

    xtr_reader = xtr.ExtractedPropertyUnstructuredGridReader()
    xtr_reader.SetInputConnection(gmy_reader.GetOutputPort())
    xtr_reader.SetExtraction(extraction)
    xtr_reader.SetTime(1000)

    xtr_reader.Update()
    out = xtr_reader.GetOutputDataObject(0)
    cd = out.GetCellData()
    assert cd.GetNumberOfArrays() == 3
    n_sites = 44250
    name_to_ncomp = {"pressure": 1, "developed_velocity_field": 3, "shearstress": 1}
    for i in range(3):
        arr = cd.GetArray(i)
        ncomp = name_to_ncomp[arr.GetName()]
        assert arr.GetNumberOfComponents() == ncomp
