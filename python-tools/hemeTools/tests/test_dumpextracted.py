from io import StringIO
import os.path

from hemeTools.converters.ExtractedPropertyTextDump import unpack


def test_simple(diffTestDir):
    # Just a smoke test
    snap = os.path.join(diffTestDir, "CleanExtracted", "flow_snapshot.dat")
    buf = StringIO()
    unpack(snap, stream=buf)
