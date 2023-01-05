# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from pathlib import Path
import numpy as np
from hlb.parsers.offset import OffsetFile


def test_load_off(diffTestDir):
    off_path = Path(diffTestDir) / "CleanExtracted" / "flow_snapshot.off"
    off = OffsetFile(off_path)
    assert off.NumberOfRanks == 2
    assert off.Data.shape == (3,)
