# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import numpy as N

class Plane(object):
    """Simple class to represent a plane.
    """
    
    def __init__(self, normal=None, centre=None, name=None, pressure=None):
        """All arguments optional. name is only used for pretty
        printing.
        
        """
        if normal is None:
            normal = (0,0,1.)
        self.normal = N.array(normal, dtype=N.float)
        if centre is None:
            centre = (0.,0.,0.)
        self.centre = N.array(centre, dtype=N.float)
        
        self.name = name
        return
    
    def __str__(self):
        if self.name is None:
            namestr = ''
        else:
            namestr = self.name + ' '
            pass
        ans  = '%snormal = %s\n' % (namestr, str(self.normal))
        ans += '%scentre = %s\n' % (namestr, str(self.centre))
        return ans
        
    def GetOrigin(self):
        return self.centre
    def GetNormal(self):
        return self.normal
    
    pass

class Pressure(object):
    """Imposed pressure on an inlet or outlet. All in mmHg.
    """
    def __init__(self, mean, amp, phase):
        """Mean pressure, amplitude of pressure, phase.
        Pressures in mmHg, phase in degrees.
        """
        self.mean = mean
        self.amp = amp
        self.phase = phase
        return
    pass
