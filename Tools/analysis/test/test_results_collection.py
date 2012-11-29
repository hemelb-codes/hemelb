# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import os
import unittest
import fixtures
from ..results_collection import ResultsCollection

class TestResultsCollection(unittest.TestCase):
	def setUp(self):
	    self.rc=ResultsCollection(fixtures.Results('cylinders').path,fixtures.ResultsConfig('example'))
	def test_construct(self):
		self.assertEqual(3*9,len(self.rc.results))
		
