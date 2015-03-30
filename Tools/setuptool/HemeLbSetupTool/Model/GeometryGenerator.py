import Generation
from .Iolets import Inlet, Outlet
from .XmlWriter import XmlWriter
from .Timer import Timer

def DVfromV(v):
    """Translate a Model.Vector.Vector to a Generation.DoubleVector.
    """
    return Generation.DoubleVector(v.x, v.y, v.z)

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
        t = Timer()
        t.Start()
        self.generator.Execute(self.skipNonIntersectingBlocks)
        XmlWriter(self).Write()
        t.Stop()
        print "Setup time: %f s" % t.GetTime()
        return

    pass
