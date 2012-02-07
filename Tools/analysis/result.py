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

class FileModel(object):
    def __init__(self,relative_path,loader):
        self.loader=loader
        self.path=relative_path
    def fullpath(self,result):
        return os.path.expanduser(os.path.join(result.path,self.path))
    def model(self,result):
        if result.files.get(self.path): return result.files.get(self.path)
        self.logger(result).debug("Loaded")
        result.files[self.path]=result.files.get(self.path) or self.loader(self.fullpath(result))
        return result.files[self.path]
    def logger(self,result):
        return logging.LoggerAdapter(logger,dict(file=self.fullpath(result)))

class ResultProperty(object):
    def __init__(self,label,memoized_file_model,parser,pattern):
        self.pattern=pattern
        self.file=memoized_file_model
        self.label=label
        self.parser=parser
    @staticmethod
    def parse_value(value):
        if value in ['None','none',None]:
            return None
        try:
            return int(value)
        except (TypeError,ValueError):
            try:
                return float(value)
            except (TypeError,ValueError):
                try:
                    return float(value.replace("_","."))
                except (TypeError,ValueError):
                    return value
    # This defines how, when an instance of this class is a property in a parent object, a value is obtained for it.
    def __get__(self,instance,owner):
        try:
            model=self.file.model(instance)
            if not model:
                raise ParseError("Bad file.")
            instance.properties[self.label]=instance.properties.get(self.label) or self.parse_value(self.parser(model,self.pattern))
            return instance.properties.get(self.label)
        except (IOError,ParseError):
            self.file.logger(instance).warning("Problem parsing value")
            return None

class ParseError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def index_parser(content,pattern):
    return content.get(pattern)
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
        logger.error
        raise ParseError("No element %s"%pattern)
    if attribute:
        try:
            return element.get(attribute)
        except AttributeError:
            raise ParseError("No attribute %s on element %s"%attribute,pattern)
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
        raise ParseError("Could not parse file.")
def nul_loader(path):
    return None
def identity_loader(path):
    return path

def result_model(config):
    class Result(object):
        """Model of a result"""
        proplist=[]
        @classmethod
        def define_file_properties(klass,config,loader,parser):
            if not config: return
            for path,data in config.iteritems():
                klass.define_properties(FileModel(path,loader),data,parser)

        @classmethod
        def define_properties(klass,file_model,data,parser):
            if not data: return
            for prop,pattern in data.iteritems():
                klass.proplist.append(prop)
                setattr(klass,prop,ResultProperty(prop,file_model,parser,pattern))
                
        def __init__(self,path):
            """path: the path to the result folder
               config: a dictionary specifying what aspects of the result folder to make into properties of the result
            """
            self.path=os.path.expanduser(path)
            self.name=os.path.basename(self.path)
            self.properties={key:None for key in self.proplist}
            self.files={}

        def datum(self,property):
            if property[0]=='(':
                # We have an expression to evaluate.
                # Do it in the binding of the current object
                try:
                    return eval(property,{},vars(self))
                except:
                    self.logger.warning("Problem handling expression %s"%property)
                    return None
            return getattr(self,property)

        def __str__(self):
            propstring=', '.join(["%s : %s"%(prop,self.datum(prop)) for prop in self.properties])
            return "Result %s: [%s]"%(self.name,propstring)


    Result.define_file_properties(config.get('yaml_files'),yaml_loader,index_parser)
    Result.define_file_properties(config.get('text_files'),text_loader,regex_parser)
    Result.define_file_properties(config.get('xml_files'),xml_loader,element_parser)
    Result.define_properties(FileModel('name',identity_loader),config.get('name_properties'),regex_parser)
    Result.define_properties(FileModel('fixed',nul_loader),config.get('fixed_properties'),identity_parser)
    return Result
