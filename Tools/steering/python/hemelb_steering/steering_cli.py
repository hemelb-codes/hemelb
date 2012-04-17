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
from config import config, config_user
from remote_hemelb import RemoteHemeLB

class Intervention(object):
    """Gather a group of results, and manipulate them in some way"""
    def __init__(self,clargs,stream=sys.stdout):
        self.stream=stream
        # By default, an unsupplied argument does not create a result property
        parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
        
        # Every property of a result is a potential command line argument argument
        for prop in config['steered_parameters']:
            parser.add_argument("--"+prop)

        # Additional possible argument, to invert the selection
        parser.add_argument("--monitor",action='store_true',default=False)
        options,extra=parser.parse_known_args(clargs)
        print(extra)
        print(vars(options))
        remote=extra[1]
        options=vars(options)
        
        print(remote)
        self.monitor=options.pop('monitor')

        config.update(config.get(remote,{}))
        config.update(config_user)
        config.update(config.get(remote,{}))

        self.port=config['port']
        self.steering_id=config['steering_id']
        self.address=config['address']
        
        # This Remote will have now connected, but not attempted to send/receive
        self.hemelb=RemoteHemeLB(port=self.port,address=self.address,steering_id=self.steering_id)
        
        for option in options:
            print("Setting %s to %s"%(option,options[option]))
            setattr(self.hemelb,option,options[option])
        
    def report(self):
        print(self.hemelb,file=self.stream)

    def act(self):        
        self.hemelb.step()
        self.report()
        if self.monitor:
            while True:
                self.hemelb.step()
                self.report()

def main():
    Intervention(sys.argv).act()

if __name__ == '__main__':
    main()
