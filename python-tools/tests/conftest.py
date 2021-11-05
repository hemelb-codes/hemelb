# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os.path

import pytest


@pytest.fixture
def diffTestDir():
    testdir = os.environ["HEMELB_TESTS_DIR"]
    diff = os.path.join(testdir, "diffTest")
    assert os.path.isdir(diff)
    return diff
