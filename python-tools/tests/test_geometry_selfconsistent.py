# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os.path
from hlb.parsers.geometry.self_consistency import CheckingLoader


def test_gmy_selfconsist(diffTestDir):
    config = os.path.join(diffTestDir, "config.xml")
    ldr = CheckingLoader(config)
    ldr.Load()

    assert not ldr.HadErrors
