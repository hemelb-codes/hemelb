#!/usr/bin/env python
# encoding: utf-8
"""
test_result.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import unittest
import fixtures
from ..result import Result

class TestResultsCollection(unittest.TestCase):
	def setUp(self):
	    self.rc=Result(fixtures.Results('cylinders').result_path(0),fixtures.ResultsConfig('example'))
	def test_construct(self):
		self.assertEqual(fixtures.Results('cylinders').result_path(0),self.rc.path)
	def test_name_property(self):
	    self.assertEqual('planck',self.rc.machine)
	def test_text_property(self):
	    self.assertEqual(1,self.rc.cores)
	    self.assertEqual('Release',self.rc.build_type)
	def test_xml_property(self):
	    self.assertEqual(2.19,self.rc.total)
	    self.assertEqual(1000,self.rc.steps)