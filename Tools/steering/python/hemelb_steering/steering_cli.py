#!/usr/bin/env python
# encoding: utf-8
"""
Command-line HemeLB steering interface.
See README for usage.
"""

from __future__ import print_function
from driver import Driver
import sys
from config import config, config_user

class Intervention(Driver):
    """Gather a group of results, and manipulate them in some way"""
    def __init__(self,clargs,stream=sys.stdout):
        Driver.__init__(self,clargs)
        self.stream=stream
        
        self.monitor=self.options.pop('monitor')
        self.show=self.options.pop('show',False)
        
        for option in self.options:
            print("Setting %s to %s"%(option,self.options[option]))
            setattr(self.hemelb,option,self.options[option])
        
    def define_args(self):
       # Every property of a result is a potential command line argument argument
       for prop in config['steered_parameters']:
           self.parser.add_argument("--"+prop)

       # Additional possible argument, to invert the selection
       self.parser.add_argument("--monitor",action='store_true',default=False)
       self.parser.add_argument("--show")
    
    def report(self):
        print(self.hemelb,file=self.stream)

    def act(self):        
        self.hemelb.step()
        if self.show:
            self.hemelb.image.pil(self.show).show()
        self.report()
        if self.monitor:
            while True:
                self.hemelb.step()
                self.report()

def main():
    Intervention(sys.argv).act()

if __name__ == '__main__':
    main()
