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
	    pass
	def test_construct(self):
		g=Graph(fixtures.GraphConfig('performance_versus_cores'))
		self.assertEqual({'type':'hemelb'},g.select)
		self.assertEqual(['cores'],g.curves)
		self.assertEqual(['total'],g.dependent)
		self.assertEqual(['steps'],g.independent)