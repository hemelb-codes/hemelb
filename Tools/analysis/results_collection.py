#!/usr/bin/env python
# encoding: utf-8
"""
results_collection.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
from result import Result

class ResultsCollection(object):
    def __init__(self,source_path,config):
        self.source_path=source_path
        # Glob over the source path results collection
        results=os.listdir(os.path.expanduser(source_path))
        self.results=[Result(os.path.join(source_path,result),config) for result in results]
    def filter(self,selection,invert=False):
        def filtration(result):
            ok=all([result.datum(prop)==value for prop,value in selection.iteritems()])
            return ok if not invert else not ok
        return filter(filtration,self.results)