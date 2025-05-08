# -*- coding: utf-8 -*-
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os.path
import pickle
from copy import copy
import numpy as np
import yaml

from vtk import vtkSTLReader

from ..Util.Observer import Observable
from .SideLengthCalculator import AverageSideLengthCalculator
from .Vector import Vector
# from .Iolets import ObservableListOfIolets, IoletLoader
from .Iolets import ObservableListOfIolets, IoletLoader, Inlet, Outlet  # ← added Inlet, Outlet

import types  # required for dynamic module creation

class FakeUnpickler(pickle.Unpickler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._HST = types.ModuleType("HemeLbSetupTool")
        self._classes = {}

    def _get_or_make_mod(self, moduleName):
        parts = moduleName.split(".")
        hst = parts.pop(0)
        assert hst == "HemeLbSetupTool"
        full = hst
        cur = self._HST
        while parts:
            name = parts.pop(0)
            if not hasattr(cur, name):
                mod = types.ModuleType(f"{full}.{name}")
                mod.__package__ = full
                full = mod.__name__
                setattr(cur, name, mod)
            cur = getattr(cur, name)
        return cur

    def _get_or_make_class(self, mod, className):
        try:
            return getattr(mod, className)
        except AttributeError:
            def __up__(this):
                return this  # no-op upgrade
            fake = type(className, (object,), {"__module__": mod.__name__, "__up__": __up__})
            setattr(mod, className, fake)
            return fake

    def find_class(self, moduleName, className):
        if moduleName.startswith("HemeLbSetupTool"):
            mod = self._get_or_make_mod(moduleName)
            return self._get_or_make_class(mod, className)
        return super().find_class(moduleName, className)



class LengthUnit(Observable):
    def __init__(self, sizeInMetres, name, abbrv):
        self.SizeInMetres = sizeInMetres
        self.Name = name
        self.Abbrv = abbrv
        return

    pass


metre = LengthUnit(1.0, "metre", "m")
millimetre = LengthUnit(1e-3, "millimetre", "mm")
micrometre = micron = LengthUnit(1e-6, "micrometre", "µm")


class Profile(Observable):
    """This class represents the parameters necessary to perform a
    setup for HemeLb and supplies the functionality to do it.

    The required parameters are below with defaults and all must be
    specified to actually create the setup files.

    """

    # Required parameters and defaults.
    _CloneOrder = ["StlFileUnitId", "StlFile", "VoxelSize"]
    _Args = {
        "StlFile": None,
        "StlFileUnitId": 1,
        "Iolets": ObservableListOfIolets(),
        "VoxelSize": 0.0,
        "TimeStepSeconds": 1e-4,
        "DurationSeconds": 5.0,
        "SeedPoint": Vector(),
        "OutputGeometryFile": None,
        "OutputXmlFile": None,
    }
    _UnitChoices = [metre, millimetre, micrometre]

    def __init__(self, **kwargs):
        """Required arguments may be set here through keyword arguments."""
        # Set attributes on this instance according to the keyword
        # args given here or the default dict if they aren't present.
        for a, default in Profile._Args.items():
            setattr(self, a, kwargs.pop(a, copy(default)))
            continue
        # Raise an error on a kwarg we don't understand
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)

        # We need a reader to get the polydata
        self.StlReader = vtkSTLReader()

        # And a way to estimate the voxel size
        self.SideLengthCalculator = AverageSideLengthCalculator()
        self.SideLengthCalculator.SetInputConnection(self.StlReader.GetOutputPort())

        # Dependencies for properties
        self.AddDependency("HaveValidStlFile", "StlFile")
        self.AddDependency("HaveValidOutputXmlFile", "OutputXmlFile")
        self.AddDependency("HaveValidOutputGeometryFile", "OutputGeometryFile")
        self.AddDependency("HaveValidSeedPoint", "SeedPoint.x")
        self.AddDependency("HaveValidSeedPoint", "SeedPoint.y")
        self.AddDependency("HaveValidSeedPoint", "SeedPoint.z")
        self.AddDependency("IsReadyToGenerate", "HaveValidStlFile")
        self.AddDependency("IsReadyToGenerate", "HaveValidOutputXmlFile")
        self.AddDependency("IsReadyToGenerate", "HaveValidOutputGeometryFile")
        self.AddDependency("IsReadyToGenerate", "HaveValidSeedPoint")
        self.AddDependency("StlFileUnit", "StlFileUnitId")
        self.AddDependency("VoxelSizeMetres", "VoxelSize")
        self.AddDependency("VoxelSizeMetres", "StlFileUnit.SizeInMetres")
        self.BoundingBoxSize = 0.0
        self.AddDependency("DefaultIoletRadius", "BoundingBoxSize")

        # When the STL changes, we should reset the voxel size and
        # update the vtkSTLReader.
        self.AddObserver("StlFile", self.OnStlFileChanged)
        return

    def UpdateAttributesBasedOnCmdLineArgs(self, cmdLineArgsDict):
        """Helper method that takes a dictionary with the arguments provided
        to the setup tool via command line (cmdLineArgsDict) and sets/updates
        the relevant class attributes.
        """
        # Some attributes need to be set in a given order to avoid side effects.
        # Set them first
        for attrName in self._CloneOrder:
            if attrName in cmdLineArgsDict:
                val = cmdLineArgsDict[attrName]
                if val is not None:
                    setattr(self, attrName, val)
                    cmdLineArgsDict.pop(attrName)

        # Set the rest
        for k, val in cmdLineArgsDict.items():
            if val is not None:
                setattr(self, k, val)

    def OnStlFileChanged(self, change):
        self.StlReader.SetFileName(self.StlFile)
        self.VoxelSize = self.SideLengthCalculator.GetOutputValue()
        self.StlReader.Update()
        surf = self.StlReader.GetOutput()
        surf.ComputeBounds()
        bounds = surf.GetBounds()
        # VTK standard bounding box
        # Compute diagonal length
        self.BoundingBoxSize = np.sqrt(
            (bounds[1] - bounds[0]) ** 2
            + (bounds[3] - bounds[2]) ** 2
            + (bounds[5] - bounds[4]) ** 2
        )
        return

    @property
    def HaveValidStlFile(self):
        """Read only property indicating if our STL file is valid."""
        return IsFileValid(self.StlFile, ext=".stl", exists=True)

    @property
    def HaveValidSeedPoint(self):
        if (
            np.isfinite(self.SeedPoint.x)
            and np.isfinite(self.SeedPoint.y)
            and np.isfinite(self.SeedPoint.z)
        ):
            return True
        return False

    @property
    def HaveValidOutputXmlFile(self):
        return IsFileValid(self.OutputXmlFile, ext=".xml")

    @property
    def HaveValidOutputGeometryFile(self):
        return IsFileValid(self.OutputGeometryFile, ext=".gmy")

    @property
    def IsReadyToGenerate(self):
        """Read only property indicating if we have enough information
        to do the setup.
        """
        if not self.HaveValidSeedPoint:
            return False
        if not self.HaveValidOutputXmlFile:
            return False
        if not self.HaveValidOutputGeometryFile:
            return False
        if not self.HaveValidStlFile:
            return False
        return True

    @property
    def StlFileUnit(self):
        return self._UnitChoices[self.StlFileUnitId]

    @property
    def VoxelSizeMetres(self):
        return self.VoxelSize * self.StlFileUnit.SizeInMetres

    @VoxelSizeMetres.setter
    def VoxelSizeMetres(self, value):
        self.VoxelSize = value / self.StlFileUnit.SizeInMetres
        return

    @property
    def DefaultIoletRadius(self):
        return self.BoundingBoxSize / 20.0

    def LoadFromFile(self, filename):
        root, ext = os.path.splitext(filename)
        if ext == ".pro":
            return self.LoadProfileV1(filename)
        elif ext == ".pr2":
            return self.LoadProfileV2(filename)
        else:
            raise ValueError("Unexpected extension on profile file: " + ext)

    def LoadProfileV2(self, filename):
        with open(filename) as f:
            state = yaml.load(f, yaml.SafeLoader)
        self._ResetPathsV2(state, filename)
        self.LoadFrom(state)
        return

    def LoadProfileV1(self, filename):
        with open(filename, "rb") as f:
            # restored = pickle.Unpickler(f, fix_imports=True).load()
            restored_fake = FakeUnpickler(f).load()
            halfway = restored_fake.__up__()  # convert fake object to a dict-like profile object
        
        # Manually assign fields instead of using CloneFrom(), to avoid constructor issues
        for attr in Profile._Args:
            val = getattr(halfway, attr, None)

            if attr == 'SeedPoint' and val is not None:
                # Upgrade legacy vector to real Vector instance
                self.SeedPoint = Vector(val.x, val.y, val.z)

            elif attr == 'Iolets' and val is not None:
                # Upgrade list of fake Inlet/Outlet objects to real ones
                real_iolets = ObservableListOfIolets()
                for io in val:
                    # Choose correct type (Inlet or Outlet)
                    inlet_or_outlet = Inlet() if getattr(io, 'Name', '').startswith("Inlet") else Outlet()

                    # Set required fields
                    inlet_or_outlet.Name = getattr(io, 'Name', None)
                    inlet_or_outlet.Radius = getattr(io, 'Radius', None)
                    inlet_or_outlet.Centre = Vector(io.Centre.x, io.Centre.y, io.Centre.z)
                    inlet_or_outlet.Normal = Vector(io.Normal.x, io.Normal.y, io.Normal.z)
                    inlet_or_outlet.Pressure = Vector(io.Pressure.x, io.Pressure.y, io.Pressure.z)

                    real_iolets.append(inlet_or_outlet)

                self.Iolets = real_iolets

            elif val is not None:
                # Default assignment for all other attributes
                setattr(self, attr, val)

        # Adjust file paths to be relative to the profile file
        self._ResetPaths(filename)

        # restored._ResetPaths(filename)
        # self.CloneFrom(restored)
        return

    def _ResetPaths(self, filename):
        # Now adjust the paths of filenames relative to the Profile file.
        # Note that this will work if an absolute path has been pickled as
        # os.path.join will discard previous path elements when it gets an
        # absolute path. (Of course, this will only work if that path is
        # correct!)
        basePath = os.path.dirname(os.path.abspath(filename))
        self.StlFile = os.path.abspath(os.path.join(basePath, self.StlFile))
        self.OutputGeometryFile = os.path.abspath(
            os.path.join(basePath, self.OutputGeometryFile)
        )
        self.OutputXmlFile = os.path.abspath(os.path.join(basePath, self.OutputXmlFile))
        return

    def _ResetPathsV2(self, state, filename):
        # Now adjust the paths of filenames relative to the Profile file.
        basePath = os.path.dirname(os.path.abspath(filename))
        state["StlFile"] = os.path.abspath(os.path.join(basePath, state["StlFile"]))
        state["OutputGeometryFile"] = os.path.abspath(
            os.path.join(basePath, state["OutputGeometryFile"])
        )
        state["OutputXmlFile"] = os.path.abspath(
            os.path.join(basePath, state["OutputXmlFile"])
        )
        return

    def Save(self, filename):
        basePath = str(os.path.dirname(filename))
        state = self.Yamlify()
        state["StlFile"] = os.path.relpath(state["StlFile"], basePath)
        state["OutputXmlFile"] = os.path.relpath(state["OutputXmlFile"], basePath)
        state["OutputGeometryFile"] = os.path.relpath(
            state["OutputGeometryFile"], basePath
        )
        with open(filename, "w") as outfile:
            yaml.dump(state, stream=outfile)

        return

    def Generate(self):
        from .OutputGeneration import PolyDataGenerator

        generator = PolyDataGenerator(self)
        generator.Execute()
        return

    def ResetVoxelSize(self, ignored=None):
        """Action to reset the voxel size to its default value."""
        self.VoxelSize = self.SideLengthCalculator.GetOutputValue()
        return

    pass


def IsFileValid(path, ext=None, exists=None):
    if not isinstance(path, str):
        return False
    if path == "":
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