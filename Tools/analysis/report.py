# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import os

class Report(object):
    def __init__(self, config, graphs):
        self.name = config['name']
        self.graphs = {
            label: graphs[label].specialise(config['defaults'], specialisation)
            for label, specialisation in config['graphs'].iteritems()
        }

    def prepare(self, results):
        for graph in self.graphs.values():
            graph.prepare(results)
    
    def write_datafiles(self):
        data_files_path = os.path.join(self.path, 'data_files')
        try:
            os.makedirs(data_files_path)
        except OSError:
            pass
        for label, graph in self.graphs.iteritems():
            out_file = open(os.path.join(data_files_path, label) + '.csv', 'w')
            graph.write_data(out_file)
            out_file.close()
    
    def write(self, report_path):
        self.path = report_path
        self.write_datafiles()
        
