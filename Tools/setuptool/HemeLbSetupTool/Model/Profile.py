import os.path
import pickle
from copy import copy
import numpy as np

from vtk import vtkSTLReader

from HemeLbSetupTool.Util.Observer import Observable, ObservableList
from HemeLbSetupTool.Model.SideLengthCalculator import AverageSideLengthCalculator
from HemeLbSetupTool.Model.Vector import Vector
from HemeLbSetupTool.Model.OutputGeneration import ConfigGenerator

import pdb
class Profile(Observable):
    """This class represents the parameters necessary to perform a
    setup for HemeLb and supplies the functionality to do it.

    The required parameters are below with defaults and all must be
    specified to actually create the setup files.
    
    """
    # Required parameters and defaults.
    _Args = {'StlFile': None,
             'Iolets': ObservableList(),
             'VoxelSize': None,
             'SeedPoint': Vector(),
             'OutputConfigFile': None,
             'OutputXmlFile': None,
             'StressType': 1}
    
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
        # And a way to estimate the voxel size
        self.sider = AverageSideLengthCalculator()
        self.sider.SetInputConnection(self.StlReader.GetOutputPort())

        # When the STL changes, we should reset the voxel size and
        # update the vtkSTLReader.
        self.AddObserver('StlFile', self.OnStlFileChanged)
        
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
        return
    
    def OnStlFileChanged(self, change):
        self.StlReader.SetFileName(self.StlFile)
        self.VoxelSize = self.sider.GetOutputValue()
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
        generator.Execute()
        return
    
    def ResetVoxelSize(self, ignored=None):
        """Action to reset the voxel size to its default value.
        """
        self.VoxelSize = self.sider.GetOutputValue()
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
