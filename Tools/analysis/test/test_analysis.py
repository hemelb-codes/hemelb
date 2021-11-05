# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os
import unittest

from .. import Analysis
import fixtures
import tempfile 

class TestAnalysis(unittest.TestCase):
    def setUp(self):
        self.config={
            'graphs': {'performance_versus_cores':fixtures.GraphConfig('performance_versus_cores')},
            'reports': {'planck_performance':fixtures.ReportConfig('planck_performance')},
            'results_path': fixtures.Results('cylinders').path,
            'reports_path': fixtures.ReportOutput('reports').path,
            'results': fixtures.ResultsConfig('example')
        }
        self.a=Analysis(self.config)
    def test_construct(self):
        self.assertEqual(Analysis,type(self.a))
    def test_prepare(self):
        self.a.load_data()
        self.a.prepare()
    def test_write(self):
        self.a.load_data()
        self.a.prepare()
        self.a.write()
