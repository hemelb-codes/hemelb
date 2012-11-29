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
