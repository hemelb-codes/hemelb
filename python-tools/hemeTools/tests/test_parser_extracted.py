import os.path

import numpy as np

from hemeTools.parsers.extraction import ExtractedProperty


def test_load_xtr(diffTestDir):
    snap_path = os.path.join(diffTestDir, "CleanExtracted", "flow_snapshot.dat")

    exp = ExtractedProperty(snap_path)

    assert exp.voxelSizeMetres == 0.0004
    assert exp.fieldCount == 3
    assert np.all(exp.times == [1000, 2000, 3000])
    N = 44250
    assert exp.siteCount == N
    data = exp.GetByTimeStep(exp.times[0])

    # ID is integer
    assert data.id.dtype == np.uint64
    assert np.all(data.id == np.arange(N))

    # Grid is (x,y,z) as 32 bit int
    assert data.grid.dtype == np.uint32
    assert data.grid.shape == (N, 3)
    assert np.all(data.grid.min(axis=0) == [1, 1, 1])
    assert np.all(data.grid.max(axis=0) == [15, 15, 250])

    # Position is (x,y,z) as C float
    assert data.position.dtype == np.float32
    assert data.position.shape == (N, 3)

    # Pressure is scalar C float
    assert data.pressure.dtype == np.float32
    assert data.pressure.shape == (N,)

    # velocity is (x,y,z) C float
    assert data["developed_velocity_field"].dtype == np.float32
    assert data.developed_velocity_field.shape == (N, 3)

    # shearstress is scalar C float
    assert data.shearstress.dtype == np.float32
    assert data.shearstress.shape == (N,)
