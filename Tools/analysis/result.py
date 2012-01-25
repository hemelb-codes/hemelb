#!/usr/bin/env python
# encoding: utf-8
"""
result.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import os
import re
import yaml
from xml.etree import ElementTree

class Result(object):
    """Model of a result"""
    def __init__(self,path,config):
        """path: the path to the result folder
           config: a dictionary specifying what aspects of the result folder to make into properties of the result
        """
        self.path=path
        self.name=os.path.basename(self.path)
        
        def index_parser(content,pattern):
            return content[pattern]
        def regex_parser(content,pattern):
            return re.search(pattern,content).groups()[0]
        def element_parser(content,pattern):
            attribute=None
            if type(pattern)==list:
                # we have a two-tuple in the yaml, the second argument is an attribute name for the element
                pattern,attribute=pattern
            element=content.find(pattern)
            if attribute:
                return element.get(attribute)
            else:
                return element.text
        def identity_parser(content,pattern):
            return pattern
            
        def yaml_loader(path):
             return yaml.load(open(path))
        def text_loader(path):
            return open(path).read()
        def xml_loader(path):
            return ElementTree.parse(path)
            
        self.define_file_properties(config.get('yaml_files'),yaml_loader,index_parser)
        self.define_file_properties(config.get('text_files'),text_loader,regex_parser)
        self.define_file_properties(config.get('xml_files'),xml_loader,element_parser)
        self.define_properties(self.name,config.get('name_properties'),regex_parser)
        self.define_properties(None,config.get('fixed_properties'),identity_parser)

    
    def define_file_properties(self,config,loader,parser):
        if not config: return
        for path,data in config.iteritems():
            content=loader(os.path.join(self.path,path))
            self.define_properties(content,data,parser)
    def define_properties(self,content,data,parser):
        if not data: return
        for prop,pattern in data.iteritems():
            setattr(self,prop,parser(content,pattern))