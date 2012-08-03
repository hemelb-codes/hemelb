#!/usr/bin/env python
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from __future__ import print_function

"""Estimate required VoxelSize, number of sites, and number of time steps, given vessel diameter"""

#Inputs
#Required
# Diameter: Vessel Size to be modelled

#Optional
# Resolution: Number of sites required across a vessel size: default 100
# Lengths: Number of lengths along the vessel to be modelled, in units of vessel diameter: default 10
# Speed: Flow speed in the vessel: default: diameter / 0.006s
# Based on: 
# Capilliary 5um : 0.79 mm/s doi 10.1016/0026-2862(81)90084-4, 
# Cerebral artery 3mm: 44 cm/s Aneths analg 1995 81; 84-9 and Stroke. 2000 Jul;31(7):1672-8
# Max mach: maximum mach number allowed, default 0.05
# Heart rate: Default 70 bpm
# Cycles: Number of heartbeats modelled: default 3

#Outputs
# VoxelSize: Voxel size needed
# Sites: count of lattice sites needed
# Steps: count of simulation steps needed

from unum.units import *
import sys

import argparse

class Estimator:
    namedSizes={
        'capillary':5.0*um,
        'mca':3.0*mm
    }
    namedSpeeds={
        'capillary':0.79*mm/s,
        'mca':44.0*cm/s
    }
    string_props=['name']
    unit_props=['rate','model_time','diameter','model_speed']
    out_props=['voxel_size','duration','min_speed_sound','step','cyclesteps','steps','sites','cores']
    float_props=['resolution','lengths','cycles','mach','speed']
    props=string_props+float_props+unit_props+out_props
    def __init__(self,name=None,diameter=None,resolution=30.0,lengths=1000.0,cycles=3.0,rate=70.0/min, 
                mach=0.1, speed=None, model_speed=10**6/s, model_time=1.0*min):
        self.resolution=resolution
        self.diameter=diameter
        self.lengths=lengths
        self.cycles=cycles
        self.rate=rate
        self.mach=mach
        self.speed=speed
        self.name=name
        self.model_speed=model_speed
        self.model_time=model_time
        if self.name==None and type(self.diameter)==type(None): #Unum raises incompatible units if compare unum to None
            self.name='mca'
        if self.name:
            self.speed=Estimator.namedSpeeds[self.name]
            self.diameter=Estimator.namedSizes[self.name]
        if type(self.speed)==type(None):
            self.speed=self.diameter/(0.006*s)
        self.voxel_size=(self.diameter/self.resolution).asUnit(m)
        self.sites=3*(self.resolution**3)*self.lengths
        self.duration=(self.cycles/self.rate).asUnit(s)
        self.min_speed_sound=self.speed/self.mach
        # speed_sound=self.VoxelSize/(time_step*sqrt(3))
        self.step=(self.voxel_size/(self.min_speed_sound*(3**(1.0/2)))).asUnit(us)
        self.steps=(self.duration/self.step)
        self.cyclesteps=(self.steps/self.cycles)
        self.cores=self.steps*self.sites/(self.model_speed*self.model_time)
        self.sites=str(self.sites/10**6)+" M"
        self.steps=str(self.steps/10**6)+" M"
    def Names(self):
        return Estimator.namedSizes.keys
    def Display(self,stream=sys.stdout):
        for prop in Estimator.props:
           print("%s : %s"%(prop,getattr(self,prop)),file=stream)
    
def main():
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    class EvalAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, eval(values))
    for prop in Estimator.unit_props:
        parser.add_argument("--"+prop,action=EvalAction)
    for prop in Estimator.float_props:
        parser.add_argument("--"+prop,type=float)
    for prop in Estimator.string_props:
        parser.add_argument("--"+prop)
    est=Estimator(**vars(parser.parse_args()))
    est.Display()
    
if __name__ == '__main__':