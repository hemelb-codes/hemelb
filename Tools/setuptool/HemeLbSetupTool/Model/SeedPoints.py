# -*- coding: utf-8 -*-

# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from copy import copy
import numpy as np

from HemeLbSetupTool.Util.Observer import Observable, ObservableListOf
from HemeLbSetupTool.Model.Vector import Vector

class AutoReg(type):
    """Quick metaclass to have subclasses register themselves on definition.
    """
    def __init__(cls, name, bases, dct):
        if not hasattr(cls, '_registry'):
            # Base - create an empty one
            cls._registry = {}
        else:
            # Derived, add to registry
            cls._registry[name] = cls
            pass
        type.__init__(cls, name, bases, dict)
        return
    
class Point(Observable):
    __metaclass__ = AutoReg
    """Represent boundary across which there can be flow.
    Do not instantiate
    """
    _Args = {'Name': 'Unknown seedpoint',
             # Initialize to the VTK defaults for now
             'Point': Vector(0., 0., 0.)}
    def __init__(self, **kwargs):
        it = self._Args.iteritems()
        for a, default in it:
            setattr(self, a,
                    kwargs.pop(a, copy(default)))
            continue
        
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)

        return

    def Yamlify(self):
        dic = Observable.Yamlify(self)
        dic['Type'] = self.__class__.__name__
        return dic
    
    pass

class SeedPoint(Point):
    pass

class ObservableListOfSeedPoints(ObservableListOf):
    ElementType = Point
    def _FindType(self, attr, source):
        cls = source['Type']
        return Point._registry[cls]
    
def SeedPointLoader(d):
    className = d.pop('Type')
    cls = Point._registry[className]
    inst = cls()
    inst.LoadFrom(d)
    return inst

