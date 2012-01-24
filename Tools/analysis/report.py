#!/usr/bin/env python
# encoding: utf-8
"""
report.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

class Report(object):
    def __init__(self,config,graphs):
        self.name=config['name']
        self.defaults=config['defaults']
        self.graphs=[graphs[label].specialise(self.defaults,specialisation) for label,specialisation in config['graphs'].iteritems()]

    def prepare(self,results):
        for graph in graphs:
            graph.prepare(results)