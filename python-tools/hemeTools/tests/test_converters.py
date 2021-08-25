import os.path

import numpy as np
from vtk import vtkPolyData, vtkUnstructuredGrid, vtkVoxel

from hemeTools.converters import ExtractedPropertyUnstructuredGridReader as xtr
from hemeTools.converters import GmyUnstructuredGridReader as gmy
from hemeTools.converters import GmyWallPointsReader as wp
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


def read_gmy_ug(path):
    reader = gmy.GmyUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutputDataObject(0)


def test_GmyUnstructuredGridReader(diffTestDir):
    check_difftest_gmy(read_gmy_ug(os.path.join(diffTestDir, "config.xml")))


def test_GmyUnstructuredGridReaderLatticeUnits(diffTestDir):
    check_difftest_gmy(read_gmy_ug(os.path.join(diffTestDir, "config.gmy")))


def read_gmy_wp(path):
    reader = wp.GmyWallPointsReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutputDataObject(0)


def check_gmy_wp(points):
    assert isinstance(points, vtkPolyData)
    n = 107888
    assert points.GetNumberOfPieces() == 1
    assert points.GetNumberOfPoints() == n
    assert points.GetNumberOfVerts() == n


def test_GmyWallPointsReader(diffTestDir):
    check_gmy_wp(read_gmy_wp(os.path.join(diffTestDir, "config.xml")))


def test_GmyWallPointsReaderLatticeUnits(diffTestDir):
    check_gmy_wp(read_gmy_wp(os.path.join(diffTestDir, "config.gmy")))


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
