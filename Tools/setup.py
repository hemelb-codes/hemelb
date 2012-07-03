#!/usr/bin/env python

from distutils.core import setup
setup(name='HemeTools',
      version='0.3',
      description='HemeLB tools',
      author='Rupert Nash',
      author_email='rupert.nash@ucl.ac.uk',
      packages=['hemeTools', 'hemeTools.converters', 'hemeTools.parsers', 'hemeTools.parsers.snapshot', 'hemeTools.parsers.geometry', 'hemeTools.surfacegenerator'],
     )
