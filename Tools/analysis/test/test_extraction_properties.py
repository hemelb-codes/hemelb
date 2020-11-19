# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os
import unittest
import fixtures

import result
Result=result.result_model(fixtures.ResultsConfig('extraction'))

class TestExtraction(unittest.TestCase):
    def setUp(self):
        self.rc=Result(fixtures.Results('extractions').result_path(0))
    def test_construct(self):
        self.assertEqual(fixtures.Results('extractions').result_path(0),self.rc.path)
    def test_name_property(self):
        self.assertEqual('cylRe1VxSz8.0e-05_80b9322046ce_localhost_4_10_10',self.rc.name)
    def test_field_count(self):
        self.assertEqual(1,self.rc.axial_field_count)
    def test_field_mean(self):
        self.assertAlmostEqual(80.0,self.rc.axial_mean_pressure,5)
    def test_extract_voxel_size(self):
        self.assertEqual(8e-05,self.rc.voxel_size)
