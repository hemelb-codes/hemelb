#!/usr/bin/env python
# encoding: utf-8
"""
test_graph.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import unittest

from ..graph import Graph
from ..results_collection import ResultsCollection

import fixtures

class TestGraph(unittest.TestCase):
    def setUp(self):
        self.g=Graph(fixtures.GraphConfig('performance_versus_cores'))
        self.g2=self.g.specialise({'select':{'machine':'planck'}},{'independent':'banana'})
        self.results=ResultsCollection(fixtures.Results('cylinders').path,fixtures.ResultsConfig('example'))
    def test_construct(self):
        self.assertEqual({'type':'hemelb'},self.g.select)
        self.assertEqual(['cores'],self.g.curves)
        self.assertEqual(['total'],self.g.dependent)
        self.assertEqual(['steps'],self.g.independent)
    def test_specialise(self):
        g3=self.g.specialise({'select':{'machine':'planck'}},{'independent':'banana'})
        # Original is unchanged
        self.assertEqual({'type':'hemelb'},self.g.select)
        self.assertEqual(['cores'],self.g.curves)
        self.assertEqual(['total'],self.g.dependent)
        self.assertEqual(['steps'],self.g.independent)
        # New one is modified
        self.assertEqual({'type':'hemelb','machine':'planck'},g3.select)
        self.assertEqual(['cores'],g3.curves)
        self.assertEqual(['total'],g3.dependent)
        self.assertEqual(['steps','banana'],g3.independent)
    def test_prepare(self):
        self.g2.prepare(self.results.results)
        print self.g2.curve_data
        self.assertEqual(27,len(self.g2.filtered_results))
        self.assertEqual(3,len(self.g2.curve_data))
        self.assertEqual(9,len(self.g2.curve_data.values()[0]))
