#!/usr/bin/env python

try:
    from Cython.Distutils import build_ext
except ImportError:
    print "Please install Cython"
    raise

from distutils.core import setup
from distutils.extension import Extension
import os

CLASSIFIERS = ["Development Status :: 5 - Production/Stable",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: BSD License",
               "Operating System :: MacOS",
               "Operating System :: POSIX :: Linux",
               "Programming Language :: Cython",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description
description = "Bessel functions of the first kind written in Cython"
fid = file('README_BESSEL.txt', 'r')
long_description = fid.read()
fid.close()
idx = max(0, long_description.find("Bessel functions of the first kind"))
long_description = long_description[idx:]

NAME                = 'jBessel'
MAINTAINER          = "Simone Manini"
MAINTAINER_EMAIL    = "simone.manini@gmail.com"
DESCRIPTION         = description
LONG_DESCRIPTION    = long_description
URL                 = "https://github.com/daron1337/jBessel"
DOWNLOAD_URL        = "http://pypi.python.org/pypi/jBessel"
LICENSE             = "BSD"
CLASSIFIERS         = CLASSIFIERS
AUTHOR              = "Simone Manini"
AUTHOR_EMAIL        = "simone.manini@gmail.com"
PLATFORMS           = "Linux/MacOsX"
ISRELEASED          = True
VERSION             = '0.1.4'

setup(name=NAME,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      classifiers=CLASSIFIERS,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      platforms=PLATFORMS,
      version=VERSION,
      cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("Bessel", ["Bessel.pyx"])]
     )