import os.path
import numpy as np

from xml.etree import ElementTree
import quantities as pq

import vtk
from vtk.util.numpy_support import vtk_to_numpy

from HemeLbSetupTool.Model.Profile import Profile
from HemeLbSetupTool.Model.Iolets import Inlet
from HemeLbSetupTool.Model.Vector import Vector

from ..utils.memo import memo_property
from ..utils.default_property import default_property
from ..utils.SortedVector import SortedVector
from ..utils.IterPairs import IterPairs

from .HemeLbParameters import HemeLbParameters, simplify
from .VesselNetwork import VesselSegment
from .PoiseuilleResistance import SpecificResistance

from ..vtkhelp.cells import IterCellPtIds

import pdb

def FindSeeds(surfaceSrc, inletPt, outletPts):
    """Given a vtkAlgorithm and a list of Inlets and Outlets, return
    the point IDs of the closest points to the inlets and outlets (as
    a vtkIdList).
    """
    surfaceSrc.Update()
    ptFinder = vtk.vtkPointLocator()
    ptFinder.SetDataSet(surfaceSrc.GetOutput())
    
    inletPtIds = vtk.vtkIdList()
    inletPtIds.InsertNextId(ptFinder.FindClosestPoint(inletPt))
    outletPtIds = vtk.vtkIdList()
    for outletPt in outletPts:
        outletPtIds.InsertNextId(ptFinder.FindClosestPoint(outletPt))
        continue
    return inletPtIds, outletPtIds

class IoletProxy(object):
    """Proxy for Model.Iolets.Iolet, delegates unknown attribute
    lookup to the proxied object.
    """
    
    def __init__(self, iolet, profile):
        self._iolet = iolet
        self._profile = profile
        
    def __getattr__(self, attr):
        # This is only called if self doesn't have the requested
        # attribute. We delegate to the real object
        return getattr(self._iolet, attr)
    
    @property
    def Centre(self):
        return self._profile.LengthUnit * Vector2Numpy(self._iolet.Centre)

    @property
    def Radius(self):
        return self._profile.LengthUnit * self._iolet.Radius

    @property
    def Normal(self):
        return Vector2Numpy(self._iolet.Normal)
    
    @Normal.setter
    def Normal(self, value):
        self._iolet.Normal = Numpy2Vector(value)
        return
    
    @memo_property
    def ClosestPointId(self):
        """Find the ID of the closest point to the iolet's centre.
        """
        return self._profile.CentreLinePointLocator.FindClosestPoint(self.Centre/self._profile.LengthUnit)

    @memo_property
    def ActualRadius(self):
        """Return the maximum inscribed sphere radius at the iolet's
        centre.
        """
        clPD = self._profile.CentreLinePolyData
        return self._profile.LengthUnit * clPD.GetPointData().GetArray('radius').GetTuple1(self.ClosestPointId)
    
    pass


class ProfileProxy(HemeLbParameters):
    """Add a few properties to the Profile. Delegates unknown
    attribute look up to the profile.
    """
    _OutputUnitsMap = HemeLbParameters._UnitQuantitiesToMap([pq.mmHg])
    
    def __init__(self, proFile):
        # Configurable parameters
        self.MaxAllowedMachNumber = 0.05
        self.InletReynoldsNumber = 300.0
        self.MinRadiusVoxels = 5.0
        
        # Read the profile
        self.FileName = proFile
        self._profile = Profile()
        self._profile.LoadFromFile(proFile)
        # Proxy its iolets
        self.Iolets = [IoletProxy(iolet, self) for iolet in self._profile.Iolets]

        # Work out where to store the centrelines
        base = os.path.splitext(self.StlFile)[0]        
        self.CentreLineFile = base + '_centrelines.vtp'

        # Sort iolets -> inlets/outlets
        inlets = []
        outlets = []
        for io in self.Iolets:
            if isinstance(io._iolet, Inlet):
                inlets.append(io)
            else:
                outlets.append(io)
                pass
            continue
        # Require only one inlet!
        assert len(inlets) == 1
        self.Inlet = inlets[0]
        # Require 1 or more outlet
        assert len(outlets) > 0
        self.Outlets = outlets
        
        return
        
    def __getattr__(self, attr):
        # This is only called if self doesn't have the requested
        # attribute. We delegate to the real Profile object
        return getattr(self._profile, attr)
    
    @property
    def LengthUnit(self):
        return self.StlFileUnit.SizeInMetres * pq.metre
    
    @memo_property
    def CentreLinePolyData(self):
        """Compute centrelines based on the profile, reusing our
        memoed copy or reading from the cache file if possible.
        """
        if (os.path.exists(self.CentreLineFile ) and 
            os.path.getmtime(self.CentreLineFile ) > os.path.getmtime(self.StlFile) and
            os.path.getmtime(self.CentreLineFile ) > os.path.getmtime(self.FileName)):
            # Cached!
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(self.CentreLineFile)
            reader.Update()
            return reader.GetOutput()

        # Have to compute it
        
        # Read the STL file
        reader = vtk.vtkSTLReader()
        reader.SetFileName(profile.StlFile)
        
        # Find the seed points for centreline calculation
        # Use points one iolet radius back along the normal.
        
        outletPts = []
        def scale(iolet):
            pt = (iolet.Centre - iolet.Radius * iolet.Normal)
            pt = pt / self.LengthUnit
            return pt.magnitude
        
        for iolet in self.Iolets:
            if isinstance(iolet._iolet, Inlet):
                inletPt = scale(iolet)
            else:
                outletPts.append(scale(iolet))
                pass
            continue
        
        srcPts, tgtPts = FindSeeds(reader, inletPt, outletPts)
        
        # Lazy import since it's so slow!
        from vmtk import vtkvmtk
        centreliner = vtkvmtk.vtkvmtkPolyDataCenterlines()
        centreliner.SetInputConnection(reader.GetOutputPort())
    
        centreliner.SetSourceSeedIds(srcPts)
        centreliner.SetTargetSeedIds(tgtPts)
        centreliner.SetRadiusArrayName("radius")
        centreliner.SetCostFunction("1/R")

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputConnection(centreliner.GetOutputPort())
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputConnection(cleaner.GetOutputPort())
        writer.SetFileName(self.CentreLineFile)
        writer.Write()
        return cleaner.GetOutput()

    @memo_property
    def CentreLinePointLocator(self):
        """A vtkAbstractPointLocator subclass that allows fast
        searching for points on the centreline.
        """
        clPtLoc = vtk.vtkOctreePointLocator()
        clPtLoc.SetDataSet(self.CentreLinePolyData)
        clPtLoc.BuildLocator()
        return clPtLoc

    @memo_property
    def CentreLinePoints(self):
        return vtk_to_numpy(self.CentreLinePolyData.GetPoints().GetData()) * self.LengthUnit

    @memo_property
    def CentreLineRadii(self):
        return vtk_to_numpy(self.CentreLinePolyData.GetPointData().GetArray('radius')) * self.LengthUnit

    @memo_property
    def VolumeEstimate(self):
        """Very rough estimate of the volume.  Assumes that the domain
        is made up of conical frustra corresponding to the line
        segments of the centrelines. Ignores any volume outside this
        (e.g. aneurysm sacs) and the overlap between successive
        segments at curves and bifurcations.
        """
        clPD = self.CentreLinePolyData
        nPts = clPD.GetNumberOfPoints()

        # List of the point IDs that have been done already.
        
        # Use this rather than the segments of the vessel tree as this
        # can cope with a network graph that has mergers as well as
        # bifurcations.
        seenPts = SortedVector(capacity=nPts, dtype=int)

        # Convert to numpy for easier access
        points = self.CentreLinePoints
        r = self.CentreLineRadii
                
        ans = 0.0 * pq.metre**3
        # Get the IDs for each path through the network
        for linePtIds in IterCellPtIds(clPD.GetLines()):
            # Iterate over adjacent IDs, each of which defines a segment.
            for i, j in IterPairs(linePtIds):
                if j in seenPts:
                    # Consider a segment done if the downstream end
                    # has been see before.
                    continue

                # Length of the segment
                delta_s = (np.sum((points[i] - points[j])**2))**0.5
                # s = distance along centreline
                # 
                #   r_i_____
                #    |      -------_____r_j
                #    |                   |
                #    |                   |
                # --s_i-----------------s_j-----> s Centreline
                # 
                #  V = \int_{s_i}^{s_j} \pi r^2 ds
                # 
                #  r(s) = (r_j - r_i) * (s - s_i) + r_i
                #         -----------
                #         (s_j - s_i)
                # 
                # so dr = (r_j - r_i) * ds
                #         -----------
                #         (s_j - s_i)
                # 
                # V = \pi (s_j - s_i) * \int_{r_i}^{r_j} r^2 dr
                #         -----------
                #         (r_j - r_i)
                #   = \pi (s_j - s_i) * (r_j^3 - r_i^3)
                #         -----------
                #         (r_j - r_i)
                delta_V = (r[i]**2 + r[i]*r[j] + r[j]**2) * np.pi * delta_s / 3.

                ans += delta_V
                # Mark the point as seen.
                seenPts.add(j)
        return ans
    
    def _AdjustIoletNormals(self):
        """Adjust the iolet normals such that they are parallel to the
        centre line at their location.
        """
        if hasattr(self, "_Adjusted") and self._Adjusted:
            return
        
        # Point IDs of them all
        closestIds = [iolet.ClosestPointId for iolet in self.Iolets]
        
        clinesPD = self.CentreLinePolyData
        points = vtk_to_numpy(clinesPD.GetPoints().GetData())

        # Going to iterate over all the paths through the system and
        # see which iolets have their closest point along each path.
        # We need to keep of which have been done to avoid repeating
        # ourselves. Then at each point we track back and forwards to
        # neighbouring points to get a better estimate of the local
        # normal. (Hence why we need the connectivity data of the
        # lines)
        adjustedPtIds = set()
        
        for linePtIds in IterCellPtIds(clinesPD.GetLines()):
            nPtsOnLine = len(linePtIds)
            # Find which iolets are closest to points on this line.
            for i, closest in enumerate(closestIds):
                match = np.where(linePtIds == closest)[0]
                if len(match) == 0:
                    # Not this iolet
                    continue
                
                # This iolet is closest to this line.
                iolet = self.Iolets[i]
                closestPointOnLine = match[0]
                
                if linePtIds[closestPointOnLine] in adjustedPtIds:
                    # But we've seen this one - next!
                    continue
                
                # Haven't seen this one, add it so don't repeat ourselves
                adjustedPtIds.add(linePtIds[closestPointOnLine])

                # Find the points connected to it on either side,
                # first as an index in this cell. (Note: make sure we
                # don't go outside the array bounds.)
                lo = max(closestPointOnLine - 1, 0)
                hi = min(closestPointOnLine + 1, nPtsOnLine - 1)
                # Then map the index-of-point-on-the-line to the point IDs
                lo = linePtIds[lo]
                hi = linePtIds[hi]
                # Now compute the normal along the line
                rLo = points[lo]
                rHi = points[hi]
                n = rHi - rLo
                n /= np.linalg.norm(n)
                # If it's the wrong way round, flip it.
                if np.dot(n, iolet.Normal) < 0:
                    n *= -1.
                    pass

                # Set in the iolet object.
                iolet.Normal = n
                # iolet.Normal.x = n[0]
                # iolet.Normal.y = n[1]
                # iolet.Normal.z = n[2]
                continue
            continue
        self._Adjusted = True
        return
        
    @property
    def InletRadius(self):
        return self.Inlet.ActualRadius
    
    @property
    @simplify
    def InletPeakVelocity(self):
        return self.KinematicViscosity * self.InletReynoldsNumber / (2. * self.InletRadius)

    @property
    @simplify
    def InletFlowRate(self):
        """In m^3 / s

        Using standard poiseuille flow results and Re = rho * U_max * dia / visc
        """
        return np.pi * self.InletReynoldsNumber * self.KinematicViscosity * self.InletRadius / 4.0
    
    @property
    def OutletPressure(self):
        return 0. * pq.mmHg

    @simplify
    def ResistanceForIds(self, idList):
        return SpecificResistance(self.CentreLinePoints[idList], self.CentreLineRadii[idList]) * self.DynamicViscosity
            
    def IndexForUnknown(self, unkn):
        return self.Unknowns.index(unkn)

    @memo_property
    def Unknowns(self):
        tmp = self.Tree
        return Unknowns

    @memo_property
    def EquationSystem(self):
        eqs = self.Tree.BuildSystem(self)
        assert len(eqs) == len(self.Unknowns)
        return eqs
    
    @memo_property
    def EquationSystemMatrix(self):
        # This needs to have uniform type, i.e. no quantities
        # Use SI base for everything.
        nUnknowns = len(self.Unknowns)
        eqs = self.EquationSystem
        ans = np.zeros((nUnknowns, nUnknowns+1))
        
        with self.SIunits():
            for i in xrange(nUnknowns):
                for j in xrange(nUnknowns):
                    ans[i,j] = self.ScaleForOutput(eqs[i].get(self.Unknowns[j], 0.))
                    continue
                ans[i, nUnknowns] = eqs[i]['rhs']
                continue
            pass
        return ans

    @memo_property
    def Matrix(self):
        return self.EquationSystemMatrix[:,:len(self.Unknowns)]
    @memo_property
    def Rhs(self):
        return self.EquationSystemMatrix[:, len(self.Unknowns)]
    
    @memo_property
    def Tree(self):
        # Before constructing the tree, must trim away the bits of
        # centreline beyond the iolets.
        inletPtId = self.Inlet.ClosestPointId
        outletPtIds = [out.ClosestPointId for out in self.Outlets]
        
        linePtIds = []
        for pathIdList in IterCellPtIds(self.CentreLinePolyData.GetLines()):
            # Find where the inlet is
            inletIndex = np.where(pathIdList == inletPtId)[0]
            assert inletIndex.shape == (1,)
            inletIndex = inletIndex[0]
            
            # Now the outlet
            # First figure out which outlet
            outletId = np.where(np.in1d(outletPtIds, pathIdList))[0]
            # Must be only one outlet.
            assert outletId.shape == (1,)            
            # Now where is this outlet's point ID in pathIdList?
            outletIndex = np.where(pathIdList == outletPtIds[outletId])[0]
            assert outletIndex.shape == (1,)
            outletIndex = outletIndex[0]
            
            linePtIds.append(pathIdList[inletIndex:outletIndex])
            continue

        # Now can split into the tree
        tree = VesselSegment(linePtIds)
        
        # To fully initialise the tree segments, need to BuildTerms
        self.Unknowns = tree.BuildTerms(flowRate=self.InletFlowRate, finalOutletPressure=self.OutletPressure)
        return tree
        
    @memo_property
    def SolutionVector(self):
        return np.linalg.solve(self.Matrix, self.Rhs)

    @simplify
    def SolutionForUnknown(self, unkn):
        i = self.IndexForUnknown(unkn)
        return self.SolutionVector[i] * unkn.units
    
    @memo_property
    def MaxVelocity(self):
        return self.Tree.MaxVelocity(self)

    @default_property
    @simplify
    def TimeStep(self):
        """MaxAllowedMachNumber = MaxVelocity / SpeedOfSound

        SpeedOfSound = VoxelSize / (sqrt(3) TimeStep)
        """
        return self.MaxAllowedMachNumber * self.VoxelSize / (np.sqrt(3) * self.MaxVelocity)

    @default_property
    def VoxelSize(self):
        return self.Tree.MinRadius(self) / self.MinRadiusVoxels
    
    @property
    def PulsatilePeriod(self):
        return 1.0 * pq.second
    
    @property
    def SimulationTime(self):
        return 10 * self.PulsatilePeriod

    @memo_property
    def MaxPathLength(self):
        return self.Tree.MaxLength(self)

    @property
    def MachNumber(self):
        return self.MaxVelocity / self.SpeedOfSound
    
    @property
    def WomersleyNumber(self):
        omega = 2.0 * np.pi / self.PulsatilePeriod
        return self.InletRadius * (omega / self.KinematicViscosity)**0.5
    
    @property
    @simplify
    def SoundTime(self):
        return self.MaxPathLength / self.SpeedOfSound
    
    @property
    @simplify
    def MomentumTime(self):
        return self.CentreLineRadii.max()**2 / self.KinematicViscosity

    @property
    @simplify
    def InletPressure(self):
        return self.SolutionForUnknown(self.Tree.InletPressure)

    @property
    def PressureDifference(self):
        return self.InletPressure - self.OutletPressure
    
    @property
    def TotalSups(self):
        return self.ScaleToLatticeUnits(self.VolumeEstimate) * self.ScaleToLatticeUnits(self.SimulationTime)

    @property
    def HectorTime(self):
        ans = pq.second * (self.TotalSups / 1e6)
        ans.units = pq.hour
        return ans

    def _UpdateProfile(self):
        self._AdjustIoletNormals()
        
        self._profile.TimeStepSeconds = self.ScaleForOutput(self.TimeStep)
        self._profile.DurationSeconds = self.ScaleForOutput(self.SimulationTime)
        
        self._profile.VoxelSize = float(self.VoxelSize / self.LengthUnit)
        return
    
    def RewriteProfile(self):
        self._UpdateProfile()
        
        base, ext = os.path.splitext(self.FileName)
        newFile = base + '_adjusted' + ext
        self.Save(newFile)
        return

    def Generate(self):
        self.RewriteProfile()
        self._profile.Generate()
        self.RewriteConfigXml()
        return
    
    def RewriteConfigXml(self):
        from HemeLbSetupTool.Model.XmlWriter import XmlWriter
        
        tree = ElementTree.parse(self.OutputXmlFile)
        root = tree.getroot()
        
        self.RewriteProperties(root)
        self.RewriteInletCondition(root)
        self.RewriteOutletConditions(root)

        XmlWriter.indent(root)
        
        tree.write(self.OutputXmlFile)
        return
    
    def RewriteInletCondition(self, root):
        condEl = root.find('inlets/inlet/condition')
        condEl.clear()
        
        condEl.set('type', 'pressure')
        condEl.set('subtype', 'cosine')
        
        inP = self.SolutionForUnknown(self.Tree.InletPressure)
        self.QuantityToXml(condEl, 'amplitude', inP)
        self.QuantityToXml(condEl, 'mean', inP)
        self.QuantityToXml(condEl, 'phase', np.pi * pq.radian)
        self.QuantityToXml(condEl, 'period', self.PulsatilePeriod)
        
        return

    def RewriteOutletConditions(self, root):
        for condEl in root.findall('outlets/outlet/condition'):
            condEl.clear()
            
            condEl.set('type', 'pressure')
            condEl.set('subtype', 'cosine')
            
            self.QuantityToXml(condEl, 'amplitude', self.OutletPressure)
            self.QuantityToXml(condEl, 'mean', self.OutletPressure)
            self.QuantityToXml(condEl, 'phase', 0.0 * pq.radian)
            self.QuantityToXml(condEl, 'period', self.PulsatilePeriod)
            continue
        return
    
    def RewriteProperties(self, root):
        # Remove any old ones
        propertiesEl = root.find('properties')
        while propertiesEl is not None:
            root.remove(propertiesEl)
            propertiesEl = root.find('properties')
            continue
        
        propertiesEl = ElementTree.SubElement(root, 'properties')

        # Get our plane centres
        planePtIds = self.Tree.SamplePtIds(50)
        nPoints = len(planePtIds)
        points = self.CentreLinePoints[planePtIds]
        # And the next point along for the normals
        nextPtIds = self.Tree.SamplePtIds(50, start=1)
        dx = self.CentreLinePoints[nextPtIds] - points
        normals = dx / np.sum(dx**2,axis=-1)[:, np.newaxis]**0.5
        # And the radius at that point, plus  a bit
        radii = self.CentreLineRadii[planePtIds]*1.1

        nDigits = int(np.ceil(np.log10(nPoints)))
        fileFmtStr = 'plane{:0%dd}.xtr' % nDigits
        periodStr = str(int(self.ScaleToLatticeUnits(self.PulsatilePeriod/25.)))
        for i in xrange(len(points)):
            poEl = ElementTree.SubElement(propertiesEl, 'propertyoutput',
                                          file=fileFmtStr.format(i),
                                          period=periodStr)
            geoEl = ElementTree.SubElement(poEl, 'geometry', type='plane')
            self.QuantityToXml(geoEl, 'point', points[i])
            self.QuantityToXml(geoEl, 'normal', normals[i])
            self.QuantityToXml(geoEl, 'radius', radii[i])

            ElementTree.SubElement(poEl, 'field', type='velocity')
            ElementTree.SubElement(poEl, 'field', type='pressure')
            
            continue
        return
    
    def QuantityToXml(self, parent, name, quantity):
        value = self.ScaleForOutput(quantity)
        if isinstance(value, np.ndarray):
            assert value.shape == (3,)
            value = '({},{},{})'.format(*value)
        else:
            value = str(value)
            pass

        try:
            units = self._OutputUnitsMap[quantity.dimensionality.simplified].symbol
        except KeyError:
            units = quantity.dimensionality.string
            pass
        
        return ElementTree.SubElement(parent, name, value=value, units=units)
    
    pass
    

def Vector2Numpy(v):
    """HemeLbSetupTool.Model.Vector.Vector -> numpy array
    """
    return np.array((v.x, v.y, v.z),dtype=float)

def Numpy2Vector(v):
    v = np.atleast_1d(v).squeeze()
    assert v.shape == (3,)
    return Vector(float(v[0]), float(v[1]), float(v[2]))

def GetPlane(iolet):
    """Return parameters for a vtkDisc + vtkTransform to represent an
    iolet, viz., centre, radius, Euler angles beta & gamma.
    """
    centre = Vec2Vec(iolet.Centre)
    nx, ny, nz = Vec2Vec(iolet.Normal)
    
    beta = np.arccos(nz)*180/np.pi
    gamma = np.arcsin(ny / np.sqrt(1-nz**2))*180/np.pi
    
    return centre, iolet.Radius, beta, gamma

if __name__ == "__main__":
    # Read command line
    import sys
    proFile = sys.argv[1]

    profile = ProfileProxy(proFile)

    ToPrint = ['PulsatilePeriod', 'MomentumTime', 'SoundTime', 'HectorTime', 'VolumeEstimate', 'PressureDifference', 'MaxVelocity', 'Tau', 'MachNumber', 'WomersleyNumber']
    def Print():
        headers = ['Name', 'Physical', 'Lattice']
        lengths = map(len, headers)
        table = []

        for attr in ToPrint:
            val = getattr(profile, attr)
            if isinstance(val, pq.Quantity):
                units = val.dimensionality.string
                if units == 'dimensionless':
                    units = ''
                phys = '{:.3g} {:s}'.format(float(val.magnitude), units)
            else:
                phys = '-'

            lat = '{:.3g}'.format(profile.ScaleToLatticeUnits(val))
            row = (attr, phys, lat)
            table.append(row)
            lengths = map(max, zip(lengths, map(len, row)))
            continue

        fmtStr = '{:%ds} | {:%ds} | {:%ds}' % tuple(lengths)
        header = fmtStr.format(*headers)
        print header
        print len(header)*'-'

        for row in table:
            print fmtStr.format(*row)

    Print()

