#!/usr/bin/env python
# encoding: utf-8
"""
analysis.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import functools

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
    
    def load_data():
        self.results=ResultsCollection(self.results_path,self.result_configuration)
        
    def prepare():
        for report in self.reports:
            report.prepare(results.results)
            
    def write():
        for report in self.reports:
            report.write(self.reports_path)

def main():
    analysis=Analysis(environment.config)
    analysis.load_data()
    analysis.prepare()
    analysis.write()

if __name__ == '__main__':
    main()

