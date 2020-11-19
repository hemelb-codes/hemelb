# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

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
