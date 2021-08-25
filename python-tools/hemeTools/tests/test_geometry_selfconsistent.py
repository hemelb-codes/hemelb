import os.path
from hemeTools.parsers.geometry.self_consistency import CheckingLoader


def test_gmy_selfconsist(diffTestDir):
    config = os.path.join(diffTestDir, "config.xml")
    ldr = CheckingLoader(config)
    ldr.Load()

    assert not ldr.HadErrors
