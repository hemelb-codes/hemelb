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

import results_collection

class TestExtraction(unittest.TestCase):
  def setUp(self):
      self.rc=results_collection.ResultsCollection(fixtures.Results('cylinders').path,fixtures.ResultsConfig('optional'))
  def test_optionally_present(self):
    for result in self.rc.results:
      if result.cores==1:
        self.assertIn('banana',result.properties)
        self.assertEqual(result.banana,1)
      else:
        self.assertNotIn('banana',result.properties)
  def test_optional_definition(self):
    for result in self.rc.results:
      if result.cores==1:
        self.assertIn('apple',result.properties)
        if result.voxelsize!=6:
          pass
          # conflicting definition case,
          # behaviour not defined
        else:
          self.assertEqual(result.apple,result.cores)
      elif result.voxelsize!=6:
        self.assertIn('apple',result.properties)
        self.assertEqual(result.apple,result.voxelsize)
      else:
        self.assertNotIn('apple',result.properties)
	      
