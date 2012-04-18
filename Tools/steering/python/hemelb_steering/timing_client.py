#!/usr/bin/env python
# encoding: utf-8
"""
curate.py

Created by James Hetherington on 2012-01-27.
Copyright (c) 2012 UCL. All rights reserved.
"""

from __future__ import print_function
from driver import Driver
from time import time
import sys

class Intervention(Driver):
    """Gather a group of results, and manipulate them in some way"""
    def __init__(self,clargs,stream=sys.stdout):
        Driver.__init__(self,clargs)
        self.stream=stream
        self.time = time()
        self.lsign=1
        self.last_step=0
        self.orbit=self.options.pop('orbit')
        
    def report(self):
        print(self.hemelb.time_step,self.hemelb.time_step-self.last_step,self.hemelb.Latitude,self.hemelb.Longitude,time()-self.time,self.hemelb.image.given_pixel_count)
        self.last_step=self.hemelb.time_step
        self.time = time()

    def act(self):
        while True:
            if self.orbit:
                self.hemelb.Latitude+=5*self.lsign
                self.hemelb.Longitude+=3
                if self.hemelb.Longitude>360:
                    self.hemelb.Longitude-=360
                if self.hemelb.Latitude>90:
                    self.hemelb.Latitude=180-self.hemelb.Latitude
                    self.hemelb.Longitude+=180
                    self.lsign=-1
                if self.hemelb.Latitude<-90:
                    self.hemelb.Latitude=-180-self.hemelb.Latitude
                    self.hemelb.Longitude+=180
                    self.lsign=+1
            self.hemelb.step()
            self.report()

    def define_args(self):
       self.parser.add_argument("--orbit",action='store_true',default=False)

def main():
    Intervention(sys.argv).act()

if __name__ == '__main__':
    main()