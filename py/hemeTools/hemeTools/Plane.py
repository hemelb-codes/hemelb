# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import numpy as np


class Plane(object):
    """Simple class to represent a plane."""

    def __init__(self, normal=None, centre=None, name=None, pressure=None):
        """All arguments optional. name is only used for pretty
        printing.

        """
        if normal is None:
            normal = (0, 0, 1.0)
        self.normal = np.array(normal, dtype=np.float)
        if centre is None:
            centre = (0.0, 0.0, 0.0)
        self.centre = np.array(centre, dtype=np.float)

        self.name = name
        return

    def __str__(self):
        if self.name is None:
            namestr = ""
        else:
            namestr = self.name + " "
            pass
        ans = "%snormal = %s\n" % (namestr, str(self.normal))
        ans += "%scentre = %s\n" % (namestr, str(self.centre))
        return ans

    def GetOrigin(self):
        return self.centre

    def GetNormal(self):
        return self.normal

    pass


class Pressure(object):
    """Imposed pressure on an inlet or outlet. All in mmHg."""

    def __init__(self, mean, amp, phase):
        """Mean pressure, amplitude of pressure, phase.
        Pressures in mmHg, phase in degrees.
        """
        self.mean = mean
        self.amp = amp
        self.phase = phase
        return

    pass
