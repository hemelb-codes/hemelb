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
        it = Iolet._Args.iteritems()
        for a, default in it:
            setattr(self, a,
                    kwargs.pop(a, copy(default)))
            continue
        
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)
        
        return

    pass

class SinusoidalPressureIolet(Iolet):
    def __init__(self, **kwargs):
        pressure = kwargs.pop('Pressure', Vector(80., 0., 0.))
        Iolet.__init__(self, **kwargs)
        self.Pressure = pressure
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
