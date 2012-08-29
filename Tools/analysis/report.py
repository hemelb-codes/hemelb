#!/usr/bin/env python
# encoding: utf-8
"""
report.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os

class Report(object):
    def __init__(self,config,graphs):
        self.name=config['name']
        self.graphs={label:graphs[label].specialise(config['defaults'],specialisation) for label,specialisation in config['graphs'].iteritems()}

    def prepare(self,results):
        for graph in self.graphs.values():
            graph.prepare(results)
    
    def write_datafiles(self):
        data_files_path=os.path.join(self.path,'data_files')
        try:
            os.makedirs(data_files_path)
        except OSError:
            pass
        for label,graph in self.graphs.iteritems():
            out_file = open(os.path.join(data_files_path,label)+'.csv','w')
            graph.write_data(out_file)
            out_file.close()
    
    def write(self,report_path):
        self.path=report_path
        self.write_datafiles()
        
