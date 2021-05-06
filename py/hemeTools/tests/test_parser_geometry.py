import os.path

import numpy as np
import yaml

from hemeTools.parsers.geometry.generic import NdIndexConverter, AllSolidBlock
from hemeTools.parsers.geometry import simple


def test_ndindexer():
    conv = NdIndexConverter((2, 3, 4))

    ijk = 0
    for idx in np.ndindex(2, 3, 4):
        assert conv.NdToOne(idx) == ijk
        assert np.all(conv.OneToNd(ijk) == idx)
        ijk += 1


def test_simple(diffTestDir):
    config_path = os.path.join(diffTestDir, "config.xml")
    loader = simple.ConfigLoader(config_path)

    gmy_path = os.path.join(diffTestDir, "config.gmy")
    assert loader.GmyFileName == gmy_path
    assert loader.Domain.VoxelSize == 4e-4

    loader.Load()

    # Overall checks
    dom = loader.Domain
    assert dom.Version == 4
    assert np.all(dom.BlockCounts == (3, 3, 32))
    n_blocks = 3 * 3 * 32
    assert dom.BlockSize == 8
    assert np.all(dom.SiteCounts == dom.BlockSize * dom.BlockCounts)

    # Block header checks
    n_sites = 44250
    assert loader.HeaderBytes == (n_blocks * 3 * 4)
    assert len(dom.BlockFluidSiteCounts) == n_blocks
    assert dom.BlockFluidSiteCounts.sum() == n_sites

    # Cylinder has radius 3mm, centred on z axis
    radius = 0.003
    # Ends specified by inlet and outlet
    pro_path = os.path.join(diffTestDir, "cyl.pr2")
    with open(pro_path) as f:
        profile = yaml.load(f, Loader=yaml.SafeLoader)
    zmin, zmax = sorted(io["Centre"]["z"] * 1e-3 for io in profile["Iolets"])

    def InsideCylinder(pt):
        x, y, z = pt
        if z < zmin or z > zmax:
            return False
        return x ** 2 + y ** 2 <= radius ** 2

    def GridToWorld(idx):
        return dom.Origin + dom.VoxelSize * idx

    for b in dom.Blocks:
        if isinstance(b, AllSolidBlock):
            bmin = GridToWorld(b.Index * dom.BlockSize)
            bmax = GridToWorld((b.Index + 1) * dom.BlockSize - 1)
            # This isn't truely correct but it will do
            assert not InsideCylinder(bmin)
            assert not InsideCylinder(bmax)
        else:
            for s in b.Sites:
                xyz = GridToWorld(s.Index)
                s.IsFluid == InsideCylinder(xyz)
