# -*- coding: utf-8 -*-
from HemeLbSetupTool.Util.Observer import Observable
from HemeLbSetupTool.Model.Vector import Vector

class Iolet(Observable):
    """Represent boundary across which there can be flow.
    """
    _Args = {'Name': None,
             'Centre': Vector(),
             'Normal': Vector(),
             'Radius': float("nan")}
    
    def __init__(self, **kwargs):
        it = Iolet._Args.iteritems()
        for a, default in it:
            setattr(self, a,
                    kwargs.pop(a, default))
            continue
        
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)
        
        return

    pass

class SinusoidalPressureIolet(Iolet):
    def __init__(self, **kwargs):
        pressure = kwargs.pop('Pressure', Vector())
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
