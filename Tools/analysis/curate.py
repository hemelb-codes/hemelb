#!/usr/bin/env python
# encoding: utf-8
"""
curate.py

Created by James Hetherington on 2012-01-27.
Copyright (c) 2012 UCL. All rights reserved.
"""
from __future__ import print_function
import sys
import os
import argparse
import shutil

from results_collection import ResultsCollection
import environment

class Curation(ResultsCollection):
    """Gather a group of results, and manipulate them in some way"""
    def __init__(self,source_path,results_config,clargs,stream=sys.stdout):
        super(Curation,self).__init__(source_path,results_config)
        parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
        for prop in self.results[0].properties:
            parser.add_argument("--"+prop)
        parser.add_argument("--invert",action='store_true',default=False)
        options,action=parser.parse_known_args(clargs)
        filtration=vars(options)
        invert=filtration.pop('invert',False)
        self.filtered_results=self.filter(filtration,invert)
        self.action=Action(action,stream)
    def act(self):
        for result in self.filtered_results:
            self.action(result)

class Action(object):
    def __init__(self,config,stream):
        self.program=config.pop(0)
        self.action=config.pop(0)
        self.arguments=config
        self.stream=stream
    def __call__(self,result):
        getattr(self,self.action)(result,*self.arguments)
    def report(self,result):
        print(result,file=self.stream)
    def name(self,result):
        print(result.name,file=self.stream)
    def accept(self,result):
        acceptance=open(os.path.join(result.path,'acceptance.txt'),'w')
        acceptance.write("Accepted: OK")
    def reject(self,result):
        acceptance=open(os.path.join(result.path,'acceptance.txt'),'w')
        acceptance.write("Accepted: NO")
    def delete(self,result):
        shutil.rmtree(result.path,ignore_errors=True)

def main():
    
	Curation(environment.config['results_path'],environment.config['results'],sys.argv).act()

if __name__ == '__main__':
	main()

