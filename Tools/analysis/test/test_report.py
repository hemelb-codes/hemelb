# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os
import unittest
import tempfile

from ..report import Report
from ..graph import Graph
from ..results_collection import ResultsCollection
import fixtures

class TestReport(unittest.TestCase):
    def setUp(self):
        self.graphs={'performance_versus_cores': Graph(fixtures.GraphConfig('performance_versus_cores'))}
        self.r=Report(fixtures.ReportConfig('planck_performance'),self.graphs)
        self.results=ResultsCollection(fixtures.Results('cylinders').path,fixtures.ResultsConfig('example'))
    def test_construct(self):
        r=self.r
        self.assertEqual('Performance on planck',r.name)
        self.assertEqual({'type':'hemelb','machine':'planck'},r.graphs['performance_versus_cores'].select)
        self.assertEqual(['cores'],r.graphs['performance_versus_cores'].curves)
        self.assertEqual(['total'],r.graphs['performance_versus_cores'].dependent)
        self.assertEqual(['sites'],r.graphs['performance_versus_cores'].independent)
    def test_prepare(self):
        self.r.prepare(self.results)
    def test_write(self):
        self.r.prepare(self.results)
        self.r.write(fixtures.ReportOutput('planck_performance').path)
        
