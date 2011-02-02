import os.path
from vtk import vtkSTLReader

from HemeLbSetupTool.Util.Observer import Observable, ObservableList
from HemeLbSetupTool.Model.SideLengthCalculator import AverageSideLengthCalculator
from HemeLbSetupTool.Model.Vector import Vector

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
             'OutputXmlFile': None}
    
    def __init__(self, **kwargs):
        """Required arguments may be set here through keyword arguments.
        """
        # Set attributes on this instance according to the keyword
        # args given here or the default dict if they aren't present.
        for a, default in Profile._Args.iteritems():
            setattr(self, a,
                    kwargs.pop(a, default))
            continue
        # Raise an error on a kwarg we don't understand
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)

        # We need a reader to get the polydata
        self.stlReader = vtkSTLReader()
        # And a way to estimate the voxel size
        self.sider = AverageSideLengthCalculator()
        self.sider.SetInputConnection(self.stlReader.GetOutputPort())

        # When the STL changes, we should reset the voxel size.
        self.AddObserver('StlFile', self.OnStlFileChanged)
        self.AddDependency('HaveValidStlFile', 'StlFile')
        return
    
    def OnStlFileChanged(self, change):
        self.stlReader.SetFileName(self.StlFile)
        self.VoxelSize = self.sider.GetOutputValue()
        return
    
    @property
    def HaveValidStlFile(self):
        """Read only property indicating if our STL file is valid.
        """
        if self.StlFile is None:
            return False
        if not os.path.exists(self.StlFile):
            return False
        if self.stlReader.GetOutput().GetNumberOfPoints() > 0:
            return True

        return False
    
    @classmethod
    def NewFromFile(cls, filename):
        new = cls()
        return new
    
    def Save(self, filename):
        return

    def Generate(self):
        return
    
    def ResetVoxelSize(self, ignored=None):
        """Action to reset the voxel size to its default value.
        """
        self.VoxelSize = self.sider.GetOutputValue()
        return
    
    pass
