# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from contextlib import contextmanager
import os.path
import pickle
import types

from ..Model.Profile import Profile
from ..Model.Vector import Vector
from ..Model.Iolets import Inlet, Outlet
from .Observer import Observable

old2new = {
    ("HemeLbSetupTool.Model.Profile", "Profile"): Profile,
    ("HemeLbSetupTool.Vector", "Vector"): Vector,
    ("HemeLbSetupTool.Iolets", "Inlet"): Inlet,
    ("HemeLbSetupTool.Iolets", "Outlet"): Outlet,
}


class FakeUnpickler(pickle.Unpickler):
    """Pickler to update old profiles.

    This Unpickler will not properly unpickle types from within the
    (removed) HemeLbSetupTools package. Instead, it will generate and
    cache modules on demand and classes within those. The types do two
    things:

    1. Hold the key/value data from the pickle

    2. Allow this to be updated to new data model via the `__up__`
    method.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._classes = {}
        self._HST = types.ModuleType("HemeLbSetupTool")

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
            real = old2new[(mod.__name__, className)]

            def __up__(this):
                o = real()
                for k, v in this.__dict__.items():
                    if hasattr(v, "__up__"):
                        v = v.__up__()
                    if isinstance(v, list):
                        v = [el.__up__() if hasattr(el, "__up__") else el for el in v]
                    setattr(o, k, v)
                return o

            fake = type(
                className,
                (object,),
                {
                    "__module__": mod.__name__,
                    "__up__": __up__,
                },
            )
            setattr(mod, className, fake)
            return fake

    def find_class(self, moduleName, className):
        if moduleName.startswith("HemeLbSetupTool"):
            mod = self._get_or_make_mod(moduleName)
            return self._get_or_make_class(mod, className)
        return super().find_class(moduleName, className)


def LoadFakeProfile(filename):
    # Fiddle an unpickler to give us a simple object with just the
    # pickled attributes set, unmodified.
    un = FakeUnpickler(open(filename, "rb"))

    fake = un.load()
    return fake


@contextmanager
def ProfileBasePath(profile, outfile):
    profile.BasePath = os.path.dirname(outfile)
    try:
        yield
    finally:
        del profile.BasePath
    return


def UpdateOutputGeometryFile(profile):
    try:
        outfile = profile.OutputConfigFile
        del profile.OutputConfigFile
        print("Info: updating from Config to Geometry")
        print('Info: old file "' + outfile + '"')
    except AttributeError:
        outfile = profile.OutputGeometryFile
        pass

    base, ext = os.path.splitext(outfile)
    if ext != ".gmy":
        outfile = base + ".gmy"
        print("Info: updating from .dat to .gmy")
        print('Info: new file "' + outfile + '"')
    profile.OutputGeometryFile = outfile
    return


def RebaseFilePath(profile, attr):
    filename = getattr(profile, attr)
    if os.path.isabs(filename):
        print(
            "Info: "
            + attr
            + " is an absolute path, truncating. Output profile will assume file is in the same directory as it."
        )
        setattr(profile, attr, os.path.basename(filename))
        pass
    return


def UpdateProfileAttributes(profile):
    UpdateOutputGeometryFile(profile)
    for pth in ["OutputGeometryFile", "OutputXmlFile", "StlFile"]:
        RebaseFilePath(profile, pth)
        continue

    # Ancient arb. params
    pulsatile_period_s = 60.0 / 70.0
    reference_pressure_mmHg = 80.0

    if hasattr(profile, "Cycles") ^ hasattr(profile, "Steps"):
        raise ValueError("Profile has only one of Cycles and Steps")

    if hasattr(profile, "Cycles") and hasattr(profile, "Steps"):
        # Really old
        # First check doesn't have new style time info
        assert not hasattr(profile, "TimeStepSeconds")
        assert not hasattr(profile, "DurationSeconds")
        profile.TimeStepSeconds = pulsatile_period_s / profile.Steps
        profile.DurationSeconds = pulsatile_period_s * profile.Cycles
        del profile.Cycles
        del profile.Steps

    unusedAttrs = list(profile.__dict__.keys())
    for attr in Profile._Args:
        try:
            unusedAttrs.remove(attr)
        except ValueError:
            print(
                "Warning: using default ("
                + str(Profile._Args[attr])
                + ') for missing attribute "'
                + attr
                + '"'
            )
            pass
        continue

    for attr in unusedAttrs:
        print('Warning: ignoring unknown attribute "' + attr + '"')
        continue

    for io in profile.Iolets:
        io.Pressure.x -= reference_pressure_mmHg
    return


def Upgrade(infilename, outfilename):
    # Load in to faked data/namespace only types
    oldProfile = LoadFakeProfile(infilename)
    # Do stuff to attributes
    UpdateProfileAttributes(oldProfile)
    # Convert to current types holding data but maybe a bit wonky.
    halfway = oldProfile.__up__()
    # Create a proper profile.
    newPro = Profile()
    newPro.CloneFrom(halfway)

    newPro.Save(outfilename)
    return


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description='Upgrade an saved profile from "config" to "geometry"'
    )
    p.add_argument("input")
    p.add_argument("output")
    args = p.parse_args()

    Upgrade(args.input, args.output)
