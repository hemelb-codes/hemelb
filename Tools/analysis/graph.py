#!/usr/bin/env python
# encoding: utf-8
"""
graph.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import copy

class Graph(object):
    def __init__(self,config):
        for prop in ['name','select','curves','dependent','independent']:
            setattr(self,prop,config[prop])
    def specialise(self,*specialisations):
        """Return a copy of this graph, with the configuration modified by some additional dicts of dicts.
        """
        result=copy.deepcopy(self)
        for specialisation in specialisations:
            if not specialisation:
                continue
            result.select.update(specialisation['select'])
            for prop in ['curves','dependent','independent']:
                if prop in specialisation:
                    getattr(result,prop).append(specialisation[prop])
        return result