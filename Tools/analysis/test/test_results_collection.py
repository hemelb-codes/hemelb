# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os
import unittest
import fixtures
from ..results_collection import ResultsCollection

class TestResultsCollection(unittest.TestCase):
	def setUp(self):
	    self.rc=ResultsCollection(fixtures.Results('cylinders').path,fixtures.ResultsConfig('example'))
	def test_construct(self):
		self.assertEqual(3*9,len(self.rc.results))
		
