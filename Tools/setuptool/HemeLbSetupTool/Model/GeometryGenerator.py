# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import Generation
from .Iolets import Inlet, Outlet
from .XmlWriter import XmlWriter
from .Timer import Timer

def DVfromV(v):
    """Translate a Model.Vector.Vector to a Generation.DoubleVector.
    """
    return Generation.DoubleVector(v.x, v.y, v.z)

def get_hwm():
    import os
    with open('/proc/{pid}/status'.format(pid=os.getpid())) as stats:
        for line in stats:
            if line.startswith('VmHWM:'):
                key, val, unit = line.split()
                return (int(val), unit)

class GeometryGenerator(object):

    def __init__(self):
        self.skipNonIntersectingBlocks = False

    def _MakeIoletProxies(self):
        # Construct the Iolet structs
        nIn = 0
        nOut = 0
        ioletProxies = []
        for io in self._profile.Iolets:
            proxy = Generation.Iolet()

            proxy.Centre = DVfromV(io.Centre) / self._profile.VoxelSize
            proxy.Normal = DVfromV(io.Normal)
            proxy.Radius = io.Radius / self._profile.VoxelSize

            if isinstance(io, Inlet):
                io.Id = proxy.Id = nIn
                proxy.IsInlet = True
                nIn += 1
            elif isinstance(io, Outlet):
                io.Id = proxy.Id = nOut
                proxy.IsInlet = False
                nOut += 1
                pass
            ioletProxies.append(proxy)
            continue
        return ioletProxies

    def _SetCommonGeneratorProperties(self):
        self.generator.SetOutputGeometryFile(
            str(self._profile.OutputGeometryFile))
        # We need to keep a reference to this to make sure it's not GC'ed
        self.ioletProxies = self._MakeIoletProxies()
        self.generator.SetIolets(self.ioletProxies)
        return

    def Execute(self):
        """Forward this to the C++ implementation.
        """
        m0 = get_hwm()
        t = Timer()
        t.Start()
        self.generator.Execute()
        XmlWriter(self).Write()
        t.Stop()
        print "Setup time: %f s" % t.GetTime()
        m1 = get_hwm()
        assert m1[1] == m0[1]
        print "Memory used by generator: {} {}".format(m1[0] - m0[0], m1[1])

        return

    pass
