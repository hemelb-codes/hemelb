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
import datetime
import functools
import subprocess
import csv
import numpy
from xml.etree import ElementTree

import logging
import environment
logger=logging.getLogger('parsing')

from helpers import *
from extraction import *

class FileModel(object):
    def __init__(self,relative_path,loader):
        self.loader=loader
        self.path=relative_path
        self.key=relative_path+loader.__name__
        
    def fullpath(self,result):
        return os.path.expanduser(os.path.join(result.path,self.path))
        
    def model(self,result):
        if result.files.get(self.key): return result.files.get(self.key)
        result.files[self.key]=result.files.get(self.key) or self.loader(self.fullpath(result))
        self.logger(result).debug("Loaded")
        return result.files[self.key]
        
    def logger(self,result):
        return logging.LoggerAdapter(logger,dict(context=self.fullpath(result)))

class ResultContent(object):
    def __init__(self,filter):
        self.filter=filter
        
    def model(self,result):
        return self.filter(result)
        
    def logger(self,result):
        return logging.LoggerAdapter(logger,dict(context=result.path))

class ResultProperty(object):
    def __init__(self,label,memoized_file_model,parser,pattern):
        self.pattern=pattern
        self.file=memoized_file_model
        self.label=label
        self.parser=parser
        
    @staticmethod
    def parse_value(value):
        if type(value)==numpy.ndarray:
          return value.tolist()
        if value in [1,0]:
          return value
        if value in ['None','none',None]:
            return None
        if value in ['True','true',True]:
            return True
        if value in ['False','false',False]:
            return False
        if type(value) in [list,float]:
            return value
        try:
            return int(str(value))
        except (TypeError,ValueError):
            try:
                return float(str(value))
            except (TypeError,ValueError):
                try:
                    return float(str(value).replace("_","."))
                except (TypeError,ValueError):
                    try:
                        return str(value).strip()
                    except AttributeError:
                        return value
                        
    def get(self,result):
        try:
            model=self.file.model(result)
            if not model:
                raise ParseError("Bad file.")
            value=self.parser(model,self.pattern)
            result.properties[self.label]=result.properties.get(self.label) or self.parse_value(value)
            return result.properties.get(self.label)
        except (IOError,ParseError, OSError) as err:
            self.file.logger(result).warning("Problem parsing value: %s"%err)
            return None
            
    # This defines how, when an instance of this class is a property in a parent object, a value is obtained for it.
    def __get__(self,instance,owner):
        return self.get(instance)

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
    
def eval_parser(content,pattern):
    try:
        # Since the properties are dynamic, they aren't in vars(self), so we have to build the binding.
        # Bind only the expressions in the pattern
        # The "content" supplied must respond to () to produce the binding
        return eval(pattern,globals(),content(pattern))
    except Exception as err:
        raise ParseError("Problem handling expression %s: %s"%(pattern,err))
        
def attribute_parser(content,pattern):
    return getattr(content,pattern)
    
def fncall_parser(content,pattern):
    out= content(pattern)
    return out
    
def column_parser(content,pattern):
    return [ResultProperty.parse_value(row[pattern]) for row in content]

def yaml_loader(path):
    file = open(path)
    result = yaml.load(file)
    file.close()
    return result
    
def text_loader(path):
    file = open(path)
    result = file.read()
    file.close()
    return result
    
def xml_loader(path):
    try:
        return ElementTree.parse(path)
    except ElementTree.ParseError:
        raise ParseError("Could not parse file.")
        
def stat_loader(path):
    return os.stat(path)


def csv_loader(path): 
    file = open(path)
    content=csv.reader(file) 
    file.close()
    return [row for row in content] 
  
def ssv_loader(path): 
    file = open(path)
    content = csv.reader(file,delimiter=' ')
    file.close()
    return [row for row in content] 

def geometry_header_loader(path):
    from hemeTools.parsers.geometry.simple import ConfigLoader
    class GeometryHeaderParsedException(BaseException):
        """Inherit from BaseException as this isn't really an error, 
        a la GeneratorExit."""
        
        pass
    
    class GeometryHeader(ConfigLoader):
        def OnEndHeader(self):
            # Abort
            raise GeometryHeaderParsedException
        def Load(self):
            try:
                ConfigLoader.Load(self)
            except GeometryHeaderParsedException:
                pass
            return
        @property
        def site_count(self):
            return sum(self.Domain.BlockFluidSiteCounts)
        @property
        def block_size(self):
            return self.Domain.BlockSize
        @property
        def block_count(self):
            return len(self.Domain.Blocks)
        pass
    gh = GeometryHeader(path)
    gh.Load()
    return gh

def null_filter(result):
    return None
    
def name_filter(result):
    return result.name
    
def binding_filter(result):
    # Return a binding suitable for use in eval, from the result
    # The object so returned must respond to () to generate the binding for an expression to be evaluated
    def binder(expression):
        terms=re.split("\W",expression)
        bindings_needed=set(terms).intersection(result.proplist)
        binding={key: getattr(result,key) for key in bindings_needed}
        binding.update(eval_helpers)
        return binding
    return binder
    
def shell_filter(result):
    return functools.partial(subprocess.check_output,cwd=os.path.expanduser(result.path))
    
def mercurial_filter(result):
    def generator(template):
        if not result.changeset: raise ParseError("No mercurial revision specified.")
        if result.changeset[-1]=='+':
            changeset=result.changeset[:-1]
        else:
            changeset=result.changeset
        try:
            return subprocess.check_output(["hg","log","-r",changeset,"--template",template],cwd=environment.localroot)
        except subprocess.CalledProcessError as err:
            raise ParseError("Problem calling mercurial: %s"%err)
    return generator

def add_properties_to_class(klass,config):
  klass.define_file_properties(config.get('yaml_files'),yaml_loader,index_parser)
  klass.define_file_properties(config.get('text_files'),text_loader,regex_parser)
  klass.define_file_properties(config.get('xml_files'),xml_loader,element_parser)
  klass.define_properties(ResultContent(name_filter),config.get('name_properties'),regex_parser)
  klass.define_properties(ResultContent(null_filter),config.get('fixed_properties'),identity_parser)
  klass.define_properties(ResultContent(binding_filter),config.get('compound_properties'),eval_parser)
  klass.define_file_properties(config.get('stat_properties'),stat_loader,attribute_parser)
  klass.define_properties(ResultContent(shell_filter),config.get('shell_properties'),fncall_parser)
  klass.define_properties(ResultContent(mercurial_filter),config.get('mercurial_properties'),fncall_parser)
  klass.define_file_properties(config.get('gmy_files'),geometry_header_loader,attribute_parser)
  klass.define_file_properties(config.get('ssv_files'),ssv_loader,column_parser)
  klass.define_file_properties(config.get('csv_files'),ssv_loader,column_parser)
  klass.define_file_properties(config.get('extraction_files'),extraction_loader,extraction_parser)

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
                if prop not in klass.proplist:
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
            
        def upgrade(self,optional_config,optional_test):
          if self.datum(optional_test):
            class new_model(self.__class__):
              pass
            add_properties_to_class(new_model,optional_config)
            self.__class__=new_model
            # empty the stored property hash to allow redefinitions
            self.properties={key:None for key in self.proplist} 
          return self

        def datum(self,property):
            """Return a property. If it is an unknown property, assume it is an anonymous compound property which wasn't stated beforehand."""
            if property in self.proplist:
                if hasattr(self, property):
                    return getattr(self,property)
                else:
                    return None
            return ResultProperty(property,ResultContent(binding_filter),eval_parser,property).get(self)

        def __str__(self):
            propstring=', '.join(["%s : %s"%(prop,self.datum(prop)) for prop in self.properties])
            return "Result %s: [%s]"%(self.name,propstring)
            
        def hash(self):
            return {prop:self.datum(prop) for prop in self.properties}
            
        def query(self,property,value):
            prop=self.datum(property)
            try:
                if value[0]=='<':
                    return prop<ResultProperty.parse_value(value[1:])
                if value[0]=='!':
                    return not prop==ResultProperty.parse_value(value[1:])
                elif value[0]=='>':
                    return prop>ResultProperty.parse_value(value[1:])
                elif value[0]=='~':
                    read_value=ResultProperty.parse_value(value[1:])
                    norm_error=abs( (read_value-prop)/(read_value+prop) )
                    return norm_error < 0.03
                else:
                    return prop==value
            except TypeError:
                return prop==value

    add_properties_to_class(Result,config)
    return Result
