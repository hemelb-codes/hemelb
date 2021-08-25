import os.path

import pytest


@pytest.fixture
def diffTestDir():
    testdir = os.environ["HEMELB_TESTS_DIR"]
    diff = os.path.join(testdir, "diffTest")
    assert os.path.isdir(diff)
    return diff
