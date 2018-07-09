#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

import numpy as np

ext_modules = [
    Extension("hemeTools.utils.xdr", ["hemeTools/utils/xdr.pyx"]),
    Extension("hemeTools.parsers.geometry.BaseSite", ["hemeTools/parsers/geometry/BaseSite.pyx"],
              include_dirs=[np.get_include(),]),
    ]

setup(name='HemeTools',
      version='0.3',
      description='HemeLB tools',
      author='Rupert Nash',
      author_email='rupert.nash@ucl.ac.uk',
      packages=['hemeTools', 'hemeTools.converters', 'hemeTools.parsers', 'hemeTools.parsers.snapshot', 'hemeTools.parsers.geometry', 'hemeTools.parsers.extraction', 'hemeTools.surfacegenerator', 'hemeTools.utils'],
      ext_modules=cythonize(ext_modules),
      zip_safe=False
     )
