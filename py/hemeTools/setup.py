#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from setuptools import setup, Extension
import Cython.Build
import numpy as np

ext_modules = [
    Extension(
        "hemeTools.utils.xdr",
        [
            "hemeTools/utils/xdr.pyx",
            "hemeTools/utils/xdr.pxd",
            "hemeTools/utils/cxdr.pxd",
        ],
    ),
    Extension(
        "hemeTools.parsers.geometry.BaseSite",
        ["hemeTools/parsers/geometry/BaseSite.pyx"],
        include_dirs=[
            np.get_include(),
        ],
    ),
]

setup(
    name="hemeTools",
    version="0.4",
    description="HemeLB tools",
    author="Rupert Nash",
    author_email="r.nash@epcc.ed.ac.uk",
    packages=[
        "hemeTools",
        "hemeTools.converters",
        "hemeTools.parsers",
        "hemeTools.parsers.snapshot",
        "hemeTools.parsers.geometry",
        "hemeTools.parsers.extraction",
        "hemeTools.surfacegenerator",
        "hemeTools.utils",
    ],
    ext_modules=ext_modules,
    cmdclass={"build_ext": Cython.Build.build_ext},
)
