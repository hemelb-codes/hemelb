# -*- coding: utf-8 -*-
from copy import copy

from HemeLbSetupTool.Util.Observer import Observable
from HemeLbSetupTool.Model.Vector import Vector

class Iolet(Observable):
    """Represent boundary across which there can be flow.
    """
    _Args = {'Name': None,
             # Initialize to the VTK defaults for now
             'Centre': Vector(0.,0.,0.),
             'Normal': Vector(0.,0.,1.),
             'Radius': 0.5}
    
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
    pass

class SinusoidalPressureIolet(Iolet):
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
