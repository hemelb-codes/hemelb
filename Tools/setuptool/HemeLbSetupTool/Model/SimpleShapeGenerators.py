from .GeometryGenerator import GeometryGenerator
from .Profile import Profile, metre
import Generation
from .Iolets import Inlet, Outlet
from .Vector import Vector

class CylinderGenerator(GeometryGenerator):

    def __init__(self, OutputGeometryFile, OutputXmlFile, VoxelSizeMetres,
                 Axis, LengthMetres, RadiusMetres,
                 InletPressure=None, OutletPressure=None):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        GeometryGenerator.__init__(self)
        self.Axis = Axis
        self.LengthMetres = LengthMetres
        self.RadiusMetres = RadiusMetres
        self.InletPressure = InletPressure
        self.OutletPressure = OutletPressure

        self._profile = Profile()
        self._profile.StlFileUnitId = Profile._UnitChoices.index(metre)
        self._profile.VoxelSize = VoxelSizeMetres
        self._profile.OutputGeometryFile = OutputGeometryFile
        self._profile.OutputXmlFile = OutputXmlFile
        self._MakeIolets()

        self.generator = Generation.CylinderGenerator()
        self._SetCommonGeneratorProperties()

        self.generator.SetCylinderLength(LengthMetres / VoxelSizeMetres)
        self.generator.SetCylinderRadius(RadiusMetres / VoxelSizeMetres)
        self.generator.SetCylinderCentre(Generation.DoubleVector(0., 0., 0.))
        self.generator.SetCylinderAxis(Generation.DoubleVector(*self.Axis))
        return

    def _MakeIolets(self):
        # Construct the Iolet structs
        inlet = Inlet()
        inlet.Centre = Vector(
            *(-0.5 * self.LengthMetres * n for n in self.Axis))
        inlet.Normal = Vector(*self.Axis)
        inlet.Radius = self.RadiusMetres
        if self.InletPressure is not None:
            inlet.Pressure = self.InletPressure
        self._profile.Iolets.append(inlet)

        outlet = Outlet()
        outlet.Centre = Vector(
            *(0.5 * self.LengthMetres * n for n in self.Axis))
        outlet.Normal = Vector(*(-n for n in self.Axis))
        outlet.Radius = self.RadiusMetres
        if self.OutletPressure is not None:
            outlet.Pressure = self.OutletPressure
        self._profile.Iolets.append(outlet)

        return

    pass

class SquareDuctGenerator(GeometryGenerator):

    def __init__(self, OutputGeometryFile, OutputXmlFile, VoxelSizeMetres,
                 OpenAxis, LengthVoxels, SideVoxels,
                 InletPressure=None, OutletPressure=None):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        GeometryGenerator.__init__(self)
        assert OpenAxis in (0, 1, 2)
        self.OpenAxis = OpenAxis
        self.LengthVoxels = LengthVoxels
        self.SideVoxels = SideVoxels
        self.Sizes = Generation.DoubleVector(SideVoxels, SideVoxels, SideVoxels)
        self.Sizes[OpenAxis] = LengthVoxels
        
        self.InletPressure = InletPressure
        self.OutletPressure = OutletPressure

        self._profile = Profile()
        self._profile.StlFileUnitId = Profile._UnitChoices.index(metre)
        self._profile.VoxelSize = VoxelSizeMetres
        self._profile.OutputGeometryFile = OutputGeometryFile
        self._profile.OutputXmlFile = OutputXmlFile
        self._MakeIolets()

        self.generator = Generation.SquareDuctGenerator()
        self._SetCommonGeneratorProperties()

        self.generator.SetOpenAxis(self.OpenAxis)
        lb = self.Sizes * -0.5
        self.generator.SetLowerBound(lb)
        ub = self.Sizes * 0.5
        self.generator.SetUpperBound(ub)
        return

    def _MakeIolets(self):
        # Construct the Iolet structs
        inlet = Inlet()
        c = [0., 0., 0.]
        c[self.OpenAxis] = -0.5 * self.LengthVoxels * self._profile.VoxelSize
        inlet.Centre = Vector(*c)
        
        n = Generation.DoubleVector()
        n[self.OpenAxis] = 1.
        inlet.Normal = Vector(n.x, n.y, n.z)
        
        inlet.Radius = self.SideVoxels * self._profile.VoxelSizeMetres
        if self.InletPressure is not None:
            inlet.Pressure = self.InletPressure
        self._profile.Iolets.append(inlet)

        outlet = Outlet()
        c = [0., 0., 0.]
        c[self.OpenAxis] = 0.5 * self.LengthVoxels * self._profile.VoxelSize
        outlet.Centre = Vector(*c)
        
        n = Generation.DoubleVector()
        n[self.OpenAxis] = -1.
        outlet.Normal = Vector(n.x, n.y, n.z)
        
        outlet.Radius = self.SideVoxels * self._profile.VoxelSizeMetres
        if self.OutletPressure is not None:
            outlet.Pressure = self.OutletPressure
        self._profile.Iolets.append(outlet)

        return

    pass
