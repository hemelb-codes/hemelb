# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from io import StringIO
import os.path

from hlb.converters.ExtractedPropertyTextDump import unpack


def test_simple(diffTestDir):
    # Just a smoke test
    snap = os.path.join(diffTestDir, "CleanExtracted", "flow_snapshot.xtr")
    buf = StringIO()
    unpack(snap, stream=buf)
