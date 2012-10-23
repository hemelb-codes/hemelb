# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import sys
import cPickle
import os.path
from contextlib import contextmanager

from HemeLbSetupTool.Model.Profile import Profile
from HemeLbSetupTool.Util.Observer import Observable

def __new__(cls, *args, **kwargs):
    obj = cls.RealClass.__new__(cls.RealClass, *args, **kwargs)
    cls.RealClass.__init__(obj)
    return obj

class ClassFinder(object):
    def __call__(self, moduleName, className):
        __import__(moduleName)
        klass = getattr(sys.modules[moduleName], className)
        
        if isinstance(klass, type) and issubclass(klass, Observable):
            return self.MakeFakeClass(klass)
    
        return klass
    
    def MakeFakeClass(self, klass):
        fake = type(klass.__name__, (object,),
                    {'RealClass': klass,
                     '__new__': __new__})
        fake.__module__ = klass.__module__
        return fake
    pass

finder = ClassFinder()
def LoadFakeProfile(filename):
    # Fiddle an unpickler to give us a simple object with just the
    # pickled attributes set, unmodified.
    un = cPickle.Unpickler(file(filename))
    un.find_global = finder
    
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
        print 'Info: updating from Config to Geometry'
        print 'Info: old file "' + outfile + '"'
    except AttributeError:
        outfile = profile.OutputGeometryFile
        pass
    
    base, ext = os.path.splitext(outfile)
    if ext != '.gmy':
        outfile = base + '.gmy'
        print 'Info: updating from .dat to .gmy'
        print 'Info: new file "' + outfile + '"'
    profile.OutputGeometryFile = outfile
    return

def RebaseFilePath(profile, attr):
    filename = getattr(profile, attr)
    if os.path.isabs(filename):
        print 'Info: ' + attr + ' is an absolute path, truncating. Output profile will assume file is in the same directory as it.'
        setattr(profile, attr, os.path.basename(filename))
        pass
    return
 
def UpdateProfileAttributes(profile):
    UpdateOutputGeometryFile(profile)
    for pth in ['OutputGeometryFile', 'OutputXmlFile', 'StlFile']:
        RebaseFilePath(profile, pth)
        continue
    
    unusedAttrs = profile.__dict__.keys()
    for attr in Profile._Args:
        try:
            unusedAttrs.remove(attr)
        except ValueError:
            print 'Warning: using default (' + str(Profile._Args[attr]) + ') for missing attribute "' + attr + '"'
            pass
        continue
    
    for attr in unusedAttrs:
        print 'Warning: ignoring unknown attribute "' + attr + '"'
        continue
    
    return

def Upgrade(infilename, outfilename):
    oldProfile = LoadFakeProfile(infilename)
    UpdateProfileAttributes(oldProfile)
    
    # Create a new, genuine profile
    newPro = Profile()
    newPro.CloneFrom(oldProfile)
    
    newPro.Save(outfilename)
    return

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description='Upgrade an saved profile from "config" to "geometry"')
    p.add_argument('input')
    p.add_argument('output')
    args = p.parse_args()

    Upgrade(args.input, args.output)
    
