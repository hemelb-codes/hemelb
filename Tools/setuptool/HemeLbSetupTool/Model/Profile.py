# -*- coding: utf-8 -*-
import os.path
import pickle
from copy import copy
import numpy as np

from vtk import vtkSTLReader, vtkTransform, vtkTransformFilter

from HemeLbSetupTool.Util.Observer import Observable, ObservableList
from HemeLbSetupTool.Model.SideLengthCalculator import AverageSideLengthCalculator
from HemeLbSetupTool.Model.Vector import Vector
from HemeLbSetupTool.Model.OutputGeneration import ConfigGenerator

#import pdb
import cProfile

class LengthUnit(Observable):
    def __init__(self, sizeInMetres, name, abbrv):
        self.SizeInMetres = sizeInMetres
        self.Name = name
        self.Abbrv = abbrv
        return 
    pass

metre = LengthUnit(1., 'metre', 'm')
millimetre = LengthUnit(1e-3, 'millimetre', 'mm')
micrometre = micron = LengthUnit(1e-6, 'micrometre', u'µm')

class Profile(Observable):
    """This class represents the parameters necessary to perform a
    setup for HemeLb and supplies the functionality to do it.

    The required parameters are below with defaults and all must be
    specified to actually create the setup files.
    
    """
    # Required parameters and defaults.
    _CloneOrder = ['StlFileUnitId']
    _Args = {'StlFile': None,
             'StlFileUnitId': 1,
             'Iolets': ObservableList(),
             'VoxelSize': 0.,
             'SeedPoint': Vector(),
             'OutputConfigFile': None,
             'OutputXmlFile': None,
             'StressType': 1}
    _UnitChoices = [metre, millimetre, micrometre]
    
    def __init__(self, **kwargs):
        """Required arguments may be set here through keyword arguments.
        """
        # Set attributes on this instance according to the keyword
        # args given here or the default dict if they aren't present.
        for a, default in Profile._Args.iteritems():
            setattr(self, a,
                    kwargs.pop(a, copy(default)))
            continue
        # Raise an error on a kwarg we don't understand
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)

        # We need a reader to get the polydata
        self.StlReader = vtkSTLReader()
        # Something to scale it to metres
        scale = self.StlFileUnit.SizeInMetres
        trans = vtkTransform()
        trans.Scale(scale, scale, scale)
        self.SurfaceSource = vtkTransformFilter()
        self.SurfaceSource.SetTransform(trans)
        self.SurfaceSource.SetInputConnection(self.StlReader.GetOutputPort())
        
        # And a way to estimate the voxel size
        self.SideLengthCalculator = AverageSideLengthCalculator()
        self.SideLengthCalculator.SetInputConnection(self.SurfaceSource.GetOutputPort())

        # Dependencies for properties
        self.AddDependency('HaveValidStlFile', 'StlFile')
        self.AddDependency('HaveValidOutputXmlFile', 'OutputXmlFile')
        self.AddDependency('HaveValidOutputConfigFile', 'OutputConfigFile')
        self.AddDependency('HaveValidSeedPoint', 'SeedPoint.x')
        self.AddDependency('HaveValidSeedPoint', 'SeedPoint.y')
        self.AddDependency('HaveValidSeedPoint', 'SeedPoint.z')
        self.AddDependency('IsReadyToGenerate', 'HaveValidStlFile')
        self.AddDependency('IsReadyToGenerate', 'HaveValidOutputXmlFile')
        self.AddDependency('IsReadyToGenerate', 'HaveValidOutputConfigFile')
        self.AddDependency('IsReadyToGenerate', 'HaveValidSeedPoint')
        self.AddDependency('StlFileUnit', 'StlFileUnitId')
        
        # When the STL changes, we should reset the voxel size and
        # update the vtkSTLReader.
        self.AddObserver('StlFile', self.OnStlFileChanged)
        
        # And when the units are changed update the transform
        self.AddObserver('StlFileUnitId', self.OnStlFileUnitIdChanged)
        return
    
    def OnStlFileChanged(self, change):
        self.StlReader.SetFileName(self.StlFile)
        self.VoxelSize = self.SideLengthCalculator.GetOutputValue()
        return
    
    def OnStlFileUnitIdChanged(self, change):
        """When we change what we think the units are in the file, this is triggered.
        The scaling of the read-in geometry is updated and the seed point and IOlet positions are suitably scaled.
        """
        # Get new length of '1' in the STL file
        scale = self.StlFileUnit.SizeInMetres
        # Get the old scaling
        oldScale = self.SurfaceSource.GetTransform().GetMatrix().GetElement(0,0)
        # Create a suitable, new Transformation
        trans = vtkTransform()
        trans.Scale(scale, scale, scale)
        
        # Scale the SeedPoint
        factor = scale/oldScale
        self.SeedPoint.x *= factor
        self.SeedPoint.y *= factor
        self.SeedPoint.z *= factor
        
        # The voxel size
        self.VoxelSize *= factor
        
        # Now any IOlet planes
        for io in self.Iolets:
            io.Centre.x *= factor
            io.Centre.y *= factor
            io.Centre.z *= factor
            io.Radius *= factor
            continue
        
        # We do this last to make sure the view reset is sensible
        # (as Modifying the SurfaceSource triggers this)
        self.SurfaceSource.SetTransform(trans)
        return
    
    @property
    def HaveValidStlFile(self):
        """Read only property indicating if our STL file is valid.
        """
        return IsFileValid(self.StlFile, ext='.stl', exists=True)

    @property
    def HaveValidSeedPoint(self):
        if np.isfinite(self.SeedPoint.x) and np.isfinite(self.SeedPoint.y) and np.isfinite(self.SeedPoint.z):
            return True
        return False
    
    @property
    def HaveValidOutputXmlFile(self):
        return IsFileValid(self.OutputXmlFile, ext='.xml')
    @property
    def HaveValidOutputConfigFile(self):
        return IsFileValid(self.OutputConfigFile, ext='.dat')
    
    @property
    def IsReadyToGenerate(self):
        """Read only property indicating if we have enough information
        to do the setup.
        """
        if not self.HaveValidSeedPoint:
            return False
        if not self.HaveValidOutputXmlFile:
            return False
        if not self.HaveValidOutputConfigFile:
            return False
        if not self.HaveValidStlFile:
            return False
        return True
    
    @property
    def StlFileUnit(self):
        return self._UnitChoices[self.StlFileUnitId]
    
    def LoadFromFile(self, filename):
        restored = pickle.Unpickler(file(filename)).load()
        self.CloneFrom(restored)
        return
    
    def Save(self, filename):
        outfile = file(filename, 'w')
        pickler = pickle.Pickler(outfile)
        pickler.dump(self)
        return
    
    def Generate(self):
        generator = ConfigGenerator(self)
#        generator.Execute()
        cProfile.runctx('generator.Execute()', globals(), locals(), 'generate.prof')
        return
    
    def ResetVoxelSize(self, ignored=None):
        """Action to reset the voxel size to its default value.
        """
        self.VoxelSize = self.SideLengthCalculator.GetOutputValue()
        return
    
    pass
    
def IsFileValid(path, ext=None, exists=None):
    if not isinstance(path, (str, unicode)):
        return False
    if path == '':
        return False
    
    if exists is not None:
        if os.path.exists(path) != exists:
            return False
        pass
    
    if ext is not None:
        ending = os.path.splitext(path)[1]
        if ending != ext:
            return False
        pass
    return True
