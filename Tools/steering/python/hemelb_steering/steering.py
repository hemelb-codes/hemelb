#!/usr/bin/env python
# encoding: utf-8
"""
Executable intended to be the beginnings of a python gui.
Not working well yet.
"""

import sys
import os
import argparse
import shutil
from config import config, config_user
from remote_hemelb import RemoteHemeLB

import Tkinter
import Image, ImageTk
from driver import Driver

class Steerer(Driver, Tkinter.Tk):

    def __init__(self, parent, clargs):
        Tkinter.Tk.__init__(self,parent)
        Driver.__init__(self,clargs)
       
        self.geometry('%dx%d' % (512, 512))
        self.update()
        self.lsign = 1
        self.act()
       
    def display(self):
        old_image = None
        pil_image = self.hemelb.image.pil('stress2')
        tkpi = ImageTk.PhotoImage(pil_image)
        label_image = Tkinter.Label(self, image=tkpi)
        label_image.place(x=0, y=0, width=pil_image.size[0], height=pil_image.size[1])
        if old_image is not None:
           old_image.destroy()
           old_image = label_image
        self.update()

    def act(self):         
        while True:
            self.hemelb.Latitude += 5 * self.lsign
            self.hemelb.Longitude += 3
            
            if self.hemelb.Longitude > 360:
                self.hemelb.Longitude -= 360
                
            if self.hemelb.Latitude > 90:
                self.hemelb.Latitude = 180 - self.hemelb.Latitude
                self.hemelb.Longitude += 180
                self.lsign = -1
                
            if self.hemelb.Latitude < -90:
                self.hemelb.Latitude = -180 - self.hemelb.Latitude
                self.hemelb.Longitude += 180
                self.lsign = +1
                
            self.hemelb.step()
            print self.hemelb.time_step
            self.display()

def main():
    Steerer(None, sys.argv).mainloop()

if __name__ == '__main__':
    main()
