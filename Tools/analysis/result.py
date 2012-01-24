#!/usr/bin/env python
# encoding: utf-8
"""
result.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import re

class Result(object):
    """Model of a result"""
    def __init__(self,path,config):
        """path: the path to the result folder
           config: a dictionary specifying what aspects of the result folder to make into properties of the result
        """
        self.path=path
        self.name=os.path.basename(self.path)
        self.xml_files=config['xml_files']
        self.text_files=config['text_files']
        self.name_properties=config['name_properties']
        print self.name
        for prop,pattern in self.name_properties.iteritems():
            setattr(self,prop,re.search(pattern,self.name).groups()[0])
        for path,data in self.text_files.iteritems():
            slurp=open(os.path.join(self.path,path)).read()
            for prop,pattern in data.iteritems():
                setattr(self,prop,re.search(pattern,slurp).groups()[0])
        