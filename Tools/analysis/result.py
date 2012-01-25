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

import logging
logger=logging.getLogger('parsing')

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
            match=re.search(pattern,content)
            if not match: return None
            return re.search(pattern,content).groups()[0]
        def element_parser(content,pattern):
            attribute=None
            if type(pattern)==list:
                # we have a two-tuple in the yaml, the second argument is an attribute name for the element
                pattern,attribute=pattern
            element=content.find(pattern)
            if element==None:
                self.logger.error("No element %s"%pattern)
                raise IOError
            if attribute:
                try:
                    return element.get(attribute)
                except AttributeError:
                    self.logger.error("No attribute %s on element %s"%attribute,pattern)
                    raise IOError
            else:
                return element.text
        def identity_parser(content,pattern):
            return pattern
            
        def yaml_loader(path):
             return yaml.load(open(path))
        def text_loader(path):
            return open(path).read()
        def xml_loader(path):
            try:
                return ElementTree.parse(path)
            except ElementTree.ParseError:
                self.logger.error("Could not parse file.")
                raise IOError
                
            
        self.define_file_properties(config.get('yaml_files'),yaml_loader,index_parser)
        self.define_file_properties(config.get('text_files'),text_loader,regex_parser)
        self.define_file_properties(config.get('xml_files'),xml_loader,element_parser)
        self.define_properties(self.name,config.get('name_properties'),regex_parser)
        self.define_properties(None,config.get('fixed_properties'),identity_parser)

    
    def define_file_properties(self,config,loader,parser):
        if not config: return
        for path,data in config.iteritems():
            fullpath=os.path.expanduser(os.path.join(self.path,path))
            self.logger=logging.LoggerAdapter(logger,dict(file=fullpath))
            try:
                content=loader(fullpath)
                self.define_properties(content,data,parser)
                self.logger.debug("Parsed OK")
            except IOError:
                self.logger.warning("Problem parsing file")
                # If the file didn't exist, or could not be parsed set the property to None.
                self.define_properties(None,data,lambda content,pattern : None)
    def define_properties(self,content,data,parser):
        if not data: return
        for prop,pattern in data.iteritems():
            setattr(self,prop,parser(content,pattern))