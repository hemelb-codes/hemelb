#!/usr/bin/env python
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

# encoding: utf-8
"""
analysis.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os

import environment
from results_collection import ResultsCollection
from graph import Graph
from report import Report

class Analysis(object):
    def __init__(self,config):
        self.results_path=config['results_path']
        self.reports_path=config['reports_path']
        self.graph_configuration=config['graphs']
        self.report_configuration=config['reports']
        self.result_configuration=config['results']
        self.graphs={label:Graph(data) for label,data in self.graph_configuration.iteritems()}
        self.reports={label:Report(data,self.graphs) for label,data in self.report_configuration.iteritems()}
    
    def load_data(self):
        self.results=ResultsCollection(self.results_path,self.result_configuration)
        
    def prepare(self):
        for report in self.reports.itervalues():
            report.prepare(self.results)
            
    def write(self):
        for label,report in self.reports.iteritems():
            report.write(os.path.expanduser(os.path.join(self.reports_path,label)))
            print "Report %s written to %s" % (report.name,report.path)

def main():
    analysis=Analysis(environment.config)
    analysis.load_data()
    analysis.prepare()
    analysis.write()

if __name__ == '__main__':
    main()
