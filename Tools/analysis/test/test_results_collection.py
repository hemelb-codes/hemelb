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
	    pass
	def test_construct(self):
		a=ResultsCollection(fixtures.Results('cylinders').path)