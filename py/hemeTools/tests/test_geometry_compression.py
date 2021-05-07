import filecmp
import os.path

import numpy as np
import yaml

from hemeTools.parsers.geometry.compression import Compressor, Decompressor


def test_decompress_recompress_same(diffTestDir, tmp_path):
    # OK take config.xml[config.gmy] --decomp--> d1.gmy --comp->
    # c1.gmy --decomp-> d2.gmy --comp-> c2.gmy
    #
    # We then require that d1 == d2 and c1 == c2
    #
    # Note we can't require c1 == orig as zlib versions etc will
    # subtly change the output.
    orig = os.path.join(diffTestDir, "config.xml")
    tmp = str(tmp_path)
    d1 = os.path.join(tmp, "d1.gmy")
    c1 = os.path.join(tmp, "c1.gmy")
    d2 = os.path.join(tmp, "d2.gmy")
    c2 = os.path.join(tmp, "c2.gmy")

    Decompressor(orig, d1).Load()
    Compressor(d1, c1).Load()
    Decompressor(c1, d2).Load()
    Compressor(d2, c2).Load()

    assert filecmp.cmp(d1, d2, shallow=False)
    assert filecmp.cmp(c1, c2, shallow=False)
