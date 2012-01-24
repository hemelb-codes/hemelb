#!/usr/bin/env python
# encoding: utf-8
"""
test_analysis.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import unittest

from ..graph import Graph

import fixtures

class TestGraph(unittest.TestCase):
    def setUp(self):
        self.g=Graph(fixtures.GraphConfig('performance_versus_cores'))
    def test_construct(self):
        self.assertEqual({'type':'hemelb'},self.g.select)
        self.assertEqual(['cores'],self.g.curves)
        self.assertEqual(['total'],self.g.dependent)
        self.assertEqual(['steps'],self.g.independent)
    def test_specialise(self):
        g2=self.g.specialise({'select':{'machine':'planck'}},{'independent':'banana'})
        # Original is unchanged
        self.assertEqual({'type':'hemelb'},self.g.select)
        self.assertEqual(['cores'],self.g.curves)
        self.assertEqual(['total'],self.g.dependent)
        self.assertEqual(['steps'],self.g.independent)
        # New one is modified
        self.assertEqual({'type':'hemelb','machine':'planck'},g2.select)
        self.assertEqual(['cores'],g2.curves)
        self.assertEqual(['total'],g2.dependent)
        self.assertEqual(['steps','banana'],g2.independent)
