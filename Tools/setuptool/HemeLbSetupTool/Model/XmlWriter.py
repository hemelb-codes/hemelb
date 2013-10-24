#
# Copyright (C) University College London, 2007-2012, all rights reserved.
#
# This file is part of HemeLB and is CONFIDENTIAL. You may not work
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#

#import numpy as np
import os.path
from xml.etree.ElementTree import Element, SubElement, ElementTree

#from .Profile import Profile, metre
from .Iolets import Inlet, Outlet
from .Vector import Vector

import pdb

class XmlWriter(object):
    VERSION = 2
    STRESSTYPE = 1
    
    def __init__(self, profile):
        self.profile = profile
        return

    @staticmethod
    def indent(elem, level=0):
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                XmlWriter.indent(elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i
        return
    
    
    def Write(self):
        root = Element('hemelbsettings', version=str(self.VERSION))
        self.DoSimulation(root)
        self.DoGeometry(root)
        self.DoIolets(root)
        self.DoVisualisation(root)
        self.DoInitialConditions(root)
        
        self.indent(root)

        xmlFile = file(self.profile.OutputXmlFile, 'wb')
        xmlFile.write('<?xml version="1.0" ?>\n')
        ElementTree(root).write(xmlFile)
        return

    def DoSimulation(self, root):
        sim = SubElement(root, 'simulation')
        QuantityElement(sim, 'step_length', self.profile.TimeStepSeconds, 's')
        QuantityElement(sim, 'steps', int(self.profile.DurationSeconds/self.profile.TimeStepSeconds), 'lattice')
        ValueElement(sim, 'stresstype', self.STRESSTYPE)
        return

    def DoGeometry(self, root):
        geom = SubElement(root, 'geometry')
        data = SubElement(geom, 'datafile')
        data.set('path', os.path.relpath(self.profile.OutputGeometryFile,
                                         os.path.split(self.profile.OutputXmlFile)[0]))
        return geom
    
    def DoInitialConditions(self, root):
        ic = SubElement(root, 'initialconditions')
        pressure = SubElement(ic, 'pressure')
        QuantityElement(pressure, 'uniform', 0.0, 'mmHg')
        return
    
    def DoProperties(self, root):
        SubElement(root, 'properties')
        return
    
    def DoIolets(self, root):
        inlets = SubElement(root, 'inlets')
        outlets = SubElement(root, 'outlets')

        for io in self.profile.Iolets:
            if isinstance(io, Inlet):
                iolet = SubElement(inlets, 'inlet')
            elif isinstance(io, Outlet):
                iolet = SubElement(outlets, 'outlet')
            else:
                continue
            
            condition = SubElement(iolet, 'condition')
            condition.text = 'You must define the condition to be enforced by hand. Sorry!'
            QuantityElement(iolet, 'normal', io.Normal, 'dimensionless')
            
            # Scale the centre to metres
            centre = Vector()
            for dim in ('x', 'y', 'z'):
                setattr(centre, dim, getattr(io.Centre, dim) * self.profile.StlFileUnit.SizeInMetres)
                
            QuantityElement(iolet, 'position', centre, 'm')
            continue
        return

    def DoVisualisation(self, root):
        vis = SubElement(root, 'visualisation')
        
        QuantityElement(vis, 'centre', Vector(0.,0.,0.), 'm')
        
        orientation = SubElement(vis, 'orientation')
        QuantityElement(orientation, 'longitude', 45.0, 'deg')
        QuantityElement(orientation, 'latitude', 45.0, 'deg')
        
        display = SubElement(vis, 'display', zoom='1.0', brightness='0.03')
        
        range = SubElement(vis, 'range')
        QuantityElement(range, 'maxvelocity', 0.1, 'm/s')
        QuantityElement(range, 'maxstress', 0.1, 'Pa')
        
        return

    pass

def ValueElement(parent, name, value):
    return SubElement(parent, name, value=str(value))

def QuantityElement(parent, name, value, units):
    return SubElement(parent, name, value=str(value), units=units)
    