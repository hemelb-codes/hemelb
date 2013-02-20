#!/usr/bin/env python
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
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
      cmdclass = {'build_ext': build_ext},
      ext_modules=ext_modules
     )
