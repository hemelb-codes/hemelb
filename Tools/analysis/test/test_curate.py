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
import StringIO
from ..curate import Curation

class TestCuration(unittest.TestCase):
    def setUp(self):
        self.buff=StringIO.StringIO()
        self.rc=Curation(fixtures.Results('cylinders').path,fixtures.ResultsConfig('example'),['myprog','name'],self.buff)
    def test_construct(self):
        self.assertEqual(27,len(self.rc.filtered_results))
        self.assertEqual(self.rc.action.action,'name')
    def test_name(self):
        print self.buff.getvalue()
        self.rc.act()
        self.assertEqual(self.buff.getvalue().split('\n')[0],'cylinder_0_001_1000_3_546058666e20_planck_1')
