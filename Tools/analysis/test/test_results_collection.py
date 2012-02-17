#!/usr/bin/env python
# encoding: utf-8
"""
test_results_collection.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import unittest
import fixtures
from ..results_collection import ResultsCollection

class TestResultsCollection(unittest.TestCase):
	def setUp(self):
	    self.rc=ResultsCollection(fixtures.Results('cylinders').path,fixtures.ResultsConfig('example'))
	def test_construct(self):
		self.assertEqual(3*9,len(self.rc.results))
		