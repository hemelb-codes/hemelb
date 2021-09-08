# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os
import sys

from skbuild import setup

if sys.platform == "darwin":
    # Python thinks it's so smart and sets the
    # MACOSX_DEPLOYMENT_TARGET environment variable that messes around
    # with what features of the compiler and C++ std lib are
    # available. Set this to use your current one.
    import platform

    release, versioninfo, machine = platform.mac_ver()
    os.environ["MACOSX_DEPLOYMENT_TARGET"] = release

setup(
    name="HemeLbSetupTool",
    version="1.2",
    author="Rupert Nash",
    author_email="r.nash@epcc.ed.ac.uk",
    packages=[
        "HemeLbSetupTool",
        "HemeLbSetupTool.Bindings",
        "HemeLbSetupTool.Util",
        "HemeLbSetupTool.Model",
        "HemeLbSetupTool.View",
        "HemeLbSetupTool.Controller",
    ],
    scripts=[
        "scripts/hemelb-setup",
        "scripts/hemelb-setup-nogui",
        "scripts/upgrade-profile",
    ],
)
