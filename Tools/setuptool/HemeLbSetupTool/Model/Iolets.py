# -*- coding: utf-8 -*-
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from copy import copy
from vtk import vtkPlaneSource
import numpy as np

from HemeLbSetupTool.Util.Observer import Observable
from HemeLbSetupTool.Model.Vector import Vector

class Iolet(Observable):
    """Represent boundary across which there can be flow.
    Do not instantiate
    """
    _Args = {'Name': 'Unknown iolet',
             # Initialize to the VTK defaults for now
             'Centre': Vector(0., 0., 0.),
             'Normal': Vector(0., 0., 1.),
             'Radius': 0.5}
    # TODO: Move the representation of the plane (a vtkPlaneSource) 
    # into this class from PlacedIolet. It should have the side 
    # effect of simplifying the bindings.
    def __init__(self, **kwargs):
        it = self._Args.iteritems()
        for a, default in it:
            setattr(self, a,
                    kwargs.pop(a, copy(default)))
            continue
        
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)

        return
    
    def __getstate__(self):
        picdic = {}
        for attr in self._Args:
            picdic[attr] = getattr(self, attr)
            continue
        return picdic
    
    @property
    def Plane(self):
        p = vtkPlaneSource()
        # Default is:
        # Centre = 0,0,0
        # Normal = 0,0,1
        # Origin = -0.5, -0.5, 0
        # So set the radius first while it's simple
        p.SetOrigin(-self.Radius, -self.Radius, 0.)
        p.SetPoint1(+self.Radius, -self.Radius, 0.)
        p.SetPoint2(-self.Radius, +self.Radius, 0.)
        # Shift to the right place
        p.SetCenter(self.Centre.x, self.Centre.y, self.Centre.z)
        # Orient correctly
        p.SetNormal(self.Normal.x, self.Normal.y, self.Normal.z)
        
        return p
    
    def IntersectWithLine(self, p1, p2):
        """Given a line defined by the two points p1,p2; and a plane defined by the
        normal n and point p0, compute an intersection. The parametric
        coordinate along the line is returned in t, and the coordinates of 
        intersection are returned in x. A zero is returned if the plane and line
        do not intersect between (0<=t<=1). If the plane and line are parallel,
        zero is returned and t is set to VTK_LARGE_DOUBLE.
        """
        tol = 1e-6
        n = np.array([self.Normal.x, self.Normal.y, self.Normal.z])
        p0 =np.array([self.Centre.x, self.Centre.y, self.Centre.z])
        # Compute line vector
        p21 = np.array(p2) - p1
    
        # Compute denominator.  If ~0, line and plane are parallel.
        num = np.dot(n,p0) - np.dot(n ,p1)
        den = np.dot(n, p21)
        
        # If denominator with respect to numerator is "zero", then the line and
        # plane are considered parallel. 
        if np.fabs(den) <= np.fabs(num*tol):
            return np.finfo(float).max, None
        
        # valid intersection
        t = num / den
        x = p1 + t * p21
        
        rSq = np.sum((x - p0)**2)
        
        if rSq <= self.Radius:
            return t, x
        
        return np.finfo(float).max, None

    pass

class SinusoidalPressureIolet(Iolet):
    """Do not instantiate
    """
    _Args = Iolet._Args.copy()
    _Args['Pressure'] = Vector(80., 0., 0.)
    
    def __init__(self, **kwargs):
        Iolet.__init__(self, **kwargs)
        
        self.AddDependency('PressureEquation', 'Pressure.x')
        self.AddDependency('PressureEquation', 'Pressure.y')
        self.AddDependency('PressureEquation', 'Pressure.z')
        return

    @property
    def PressureEquation(self):
        try:
            avg = self.Pressure.x
            amp = self.Pressure.y
            phs = self.Pressure.z
            ans = u'p = %.2f + %.2f cos(wt + %.0f°)' % (avg, amp, phs)
            return ans
        except:
            return ''
        return
    
    pass

class Inlet(SinusoidalPressureIolet):
    pass

class Outlet(SinusoidalPressureIolet):
    pass
