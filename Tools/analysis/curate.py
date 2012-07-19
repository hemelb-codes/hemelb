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
import csv
import cProfile
from pprint import PrettyPrinter

from results_collection import ResultsCollection
from result import ResultProperty
import environment

class Curation(ResultsCollection):
    """Gather a group of results, and manipulate them in some way"""
    def __init__(self,source_path,results_config,clargs,stream=sys.stdout):
        # Prepare to parse the arguments
        super(Curation,self).__init__(source_path,results_config)
        # By default, an unsupplied argument does not create a result property
        parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)

        # The parser uses Result.parse_value to generate float or int values.
        class ParseAction(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None):
                setattr(namespace, self.dest, ResultProperty.parse_value(values))

        # Every property of a result is a potential command line argument argument
        for prop in self.results[0].properties:
            parser.add_argument("--"+prop,action=ParseAction)

        # Additional possible argument, to invert the selection
        parser.add_argument("--invert",action='store_true',default=False)
        parser.add_argument("--latest",type=int)
        options,action=parser.parse_known_args(clargs)

        # Ok, we have the results of argument parsing, get a dictionary of them and use it to filter the results
        filtration=vars(options)
        filter_opts={}
        filter_opts['invert']=filtration.pop('invert',False)
        filter_opts['latest']=filtration.pop('latest',False)
        self.filtered_results=self.filter(filtration,**filter_opts)
        self.action=Action(action,stream)

    def act(self):
        self.action.start()
        for result in self.filtered_results:
            self.action(result)

class Action(object):
    def __init__(self,config,stream):
        self.program=config.pop(0)
        self.action=config.pop(0)
        self.arguments=config
        self.stream=stream
        self.writer=csv.writer(self.stream, delimiter=' ')
        self.pp=PrettyPrinter(stream=self.stream)
    def start(self):
        if self.action in ['display','zip']:
            print("#",end='',file=self.stream)
            print(*self.arguments,file=self.stream)
            print("",file=self.stream)
    def __call__(self,result):
        getattr(self,self.action)(result,*self.arguments)
    def report(self,result):
        print(result,file=self.stream)
    def inspect(self,result):
        self.pp.pprint(result.hash())
    def display(self,result,*cols):
        self.writer.writerow(map(lambda x: x if not x==None else 'None',[result.datum(col) for col in cols]))
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
    def cat(self,result,*files):
        for afile in files:
            content=open(os.path.join(result.path,afile)).read()
            print(content,file=self.stream)
    def zip(self,result,*cols):
        print("#%s"%result.name, file=self.stream)
        self.writer.writerows(
            zip(*[result.datum(col) for col in cols])
        )
        print("\n",file=self.stream)

def main():
    Curation(environment.config['results_path'],environment.config['results'],sys.argv).act()

if __name__ == '__main__':
    main()
    #cProfile.run('main()','curate.prof')
