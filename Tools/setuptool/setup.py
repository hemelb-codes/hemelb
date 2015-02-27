# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import os.path
import sys
import platform
import re

import numpy
import vtk
try:
    import wxversion
    wxversion.select('2.8')
except wxversion.VersionError:
    wxversion.select('2.9')
except ImportError:
    pass

from distutils.core import setup
from distutils.extension import Extension

# from Cython.Distutils import build_ext
# from Cython.Distutils.extension import Extension

def GetVtkVersion():
    v = vtk.vtkVersion()
    major = v.GetVTKMajorVersion()
    minor = v.GetVTKMinorVersion()
    return major, minor

def GetVtkLibDirDarwin(aVtkSharedLibrary):
    def IterRpaths(obj):
        """Iterator that yields the @rpaths added to the object.
        """
        # Run otool -l $obj.

        # Find each load command and check if it's
        # 'cmd LC_RPATH'. If so, it should be followed by
        # 'cmdsize <integer>' and then by 'path <path> (offset
        # <integer>)'. Note there may be zero or more of these!
        cmd = "otool -l %s" % obj
        lines = os.popen(cmd).read().split('\n')

        # Very simple state machine to parse the output of otool
        state = 'top'
        for line in lines:
            if state == 'top':
                if line.startswith('Load command'):
                    state = 'start command'
                    pass
                continue

            if state == 'start command':
                cmd, cmdName = line.split()
                assert cmd == 'cmd', "error parsing otool output, expected 'cmd $CMD_NAME'"
                if cmdName == 'LC_RPATH':
                    state = 'in rpath size'
                else:
                    state = 'in other'
                    pass
                continue

            if state == 'in other':
                if line.startswith('Load command'):
                    state = 'start command'
                    pass
                continue

            if state == 'in rpath size':
                state = 'in rpath path'
                continue

            if state == 'in rpath path':
                cmd, path, offset = line.split(None, 2)
                assert cmd == 'path'
                yield path
                state = top
                continue
        return
    
    def DealWithLinkOptions(libDir, libBase, loaderPaths):
        
        if libDir.startswith('/'):
            # It's an absolute path - easy!
            if os.path.exists(os.path.join(libDir, libBase)):
                return os.path.normpath(libDir)
            raise ValueEror("Can't find library")
            
        elif libDir.startswith('@loader_path'):
            # Path relative to the library or executable
            for pth in loaderPaths:
                ans = libDir.replace('@loader_path', pth)
                if os.path.exists(os.path.join(ans,libBase)):
                    return os.path.normpath(ans)
                continue
            raise ValueError("Can't find library '%s'" % libDir)
        
        elif libDir.startswith('@rpath'):
            # Need to examine the rpaths set in the executable and
            # shared libraries. Potentially all the loaded ones, but
            # python and the vtk extension *should* be enough.
            for obj in (sys.executable, aVtkSharedLibrary):
                for rpath in IterRpaths(obj):
                    # replace @rpath with <path> and try to search again
                    try:
                        return DealWithLinkOptions(
                            libDir.replace("@rpath", rpath),
                            libBase, loaderPaths + (os.path.dirname(sys.executable),)
                        )
                    except ValueError:
                        pass
                    continue
                continue
            raise ValueError("Could not find actual directory in RPATH for '%s'" % libDir)
        elif libDir == '':
            # On DYLD_LIBRARY_PATH
            dyld_library_path = os.environ.get('DYLD_LIBRARY_PATH', '')
            for pth in dyld_library_path.split(':'):
                if pth == '': continue
                potential = os.path.join(pth, libBase)
                if os.path.exists(potential):
                    return os.path.normpath(pth)
                continue
            # Wasn't there....
            pass
        # Don't know how to deal with this one
        raise ValueError("Can't deduce actual directory from otool output '{}'".format(line))

    # End of helper functions.
    results = os.popen('otool -L %s' % aVtkSharedLibrary).read()
    lines = results.split('\n')
    # Ditch first line
    extLibPath = lines.pop(0)
    extLibDir = os.path.dirname(extLibPath)

    for line in lines:
        # Split on whitespace
        words = line.split()
        try:
            # First word is the full path of the library
            libPath = words[0]
        except IndexError:
            # Must have been empty
            continue
        
        try:
            # pdb.set_trace()
            libDir, libBase = os.path.split(libPath)
            if libBase.startswith('libvtkCommon'):
                return DealWithLinkOptions(libDir, libBase, (extLibDir,))
            
        except ValueError:
            pass
        continue
    # Searched them all without finding a match!
    raise ValueError("Can't figure out the VTK include dir")

def GetVtkLibDirLinux(aVtkSharedLibrary):
    results = os.popen('ldd %s' % aVtkSharedLibrary).read()

    lines = results.split('\n')
    
    for line in lines:
        try:
            # Get the second half, showing where the symbol found
            name, pathAndStuff = line.split('=>')
            
            # Split on whitespace
            words = pathAndStuff.split()
            # First word is the full path of the library
            libPath = words[0]
            dir, base = os.path.split(libPath)
            # Do a string split rather than os.path.splitext as the latter splits on the last dot
            base, rest = base.split('.', 1)
            if base == 'libvtkCommon':
                return dir
            
        except ValueError:
            pass
        continue
    # Searched them all without finding a match!
    raise ValueError("Can't figure out the VTK include dir")

def LibToInclude(vtkLibDir):
    libDir, vtk = os.path.split(vtkLibDir)
    if vtk.startswith('vtk'):
        prefix, lib = os.path.split(libDir)
    elif vtk == 'lib':
        prefix = libDir
        libDir = vtkLibDir
    else:
        raise ValueError("Can't deduce path to VTK include dir from library dir '%s'" % vtkLibDir)

    includeDir = os.path.join(prefix, 'include')
    vtkIncludeDir = os.path.join(includeDir, 'vtk-%d.%d' % GetVtkVersion())
    return vtkIncludeDir

def GetVtkLibDir():
    try:
        aVtkSharedLibrary = vtk.libvtkCommonPython.__file__
    except:
        try:
            aVtkSharedLibrary = vtk.vtkCommonPython.__file__
        except:
            try:
                aVtkSharedLibrary = vtk.vtkCommonCorePython.__file__ # Separated in the VTK "modularisation"
            except:
                import vtkCommonPython
                aVtkSharedLibrary = vtkCommonPython.__file__
    osName = platform.system()
    if osName == 'Darwin':
        return GetVtkLibDirDarwin(aVtkSharedLibrary)
    elif osName == 'Linux':
        return GetVtkLibDirLinux(aVtkSharedLibrary)
    else:
        raise ValueError("Don't know how to determing VTK path on OS '%s'" % osName)

def GetBoostDir(hemeLbDir):
    boostDir = os.path.join(hemeLbDir, '../dependencies/include/') 
    return os.path.normpath(boostDir)

def GetVtkCompileFlags(vtkLibDir):
    # SET(VTK_REQUIRED_CXX_FLAGS " -Wno-deprecated -no-cpp-precomp")
    flagFinder = re.compile(r'SET\(VTK_REQUIRED_CXX_FLAGS "(.*)"\)')
    
    try:
        cmake = file(os.path.join(vtkLibDir, 'VTKConfig.cmake'))
    except:
        return []
    
    for line in cmake:
        match = flagFinder.search(line)
        if match:
            flags = match.group(1)
            return flags.split()
        continue
    
    return []

def HaveXdrUint():
    import tempfile
    import shutil
    import subprocess
    prog = '''
#include <stdint.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
int main(int count, char** v){
  char buffer[15] = \"aaaaaaaaaaaaa\";
  XDR xdr;
  xdrmem_create(&xdr, buffer, 32, XDR_ENCODE);
  uint16_t a;
  uint32_t b;
  uint64_t c;
  xdr_uint16_t(&xdr, &a);
  xdr_uint32_t(&xdr, &b);
  xdr_uint64_t(&xdr, &c);
  return b;
}
'''
    tDir = tempfile.mkdtemp()
    try:
        curDir = os.getcwd()
        os.chdir(tDir)
        try:
            with file('test.cpp', 'w') as cSrc:
                cSrc.write(prog)
            compile = 'g++ -c test.cpp 2> /dev/null > /dev/null'
            returnCode = subprocess.call(compile, shell=True)
        finally:
            os.chdir(curDir)
    finally:
        shutil.rmtree(tDir)
    return (returnCode == 0)
    
def GetHemeLbCompileFlags():
    osName = platform.system()
    flags = []
    if osName == 'Darwin':
        flags.append('-DHEMELB_CFG_ON_BSD')
    if not HaveXdrUint():
        flags.append('-Dxdr_uint16_t=xdr_u_int16_t')
        flags.append('-Dxdr_uint32_t=xdr_u_int32_t')
        flags.append('-Dxdr_uint64_t=xdr_u_int64_t')
    return flags
    
if __name__ == "__main__":
    # numpy, vtk
    HemeLbDir = os.path.abspath('../../Code')
    BoostDir = GetBoostDir(HemeLbDir)
    vtkLibDir = GetVtkLibDir()
    if os.getenv('VTKINCLUDE'):
        vtkIncludeDir = os.getenv('VTKINCLUDE')
    else:
        vtkIncludeDir = LibToInclude(vtkLibDir)
    include_dirs = [vtkIncludeDir, HemeLbDir, BoostDir]
    
    libraries = []
    library_dirs = []
    extra_compile_args = GetVtkCompileFlags(vtkLibDir) + GetHemeLbCompileFlags()
    extra_link_args = ['-lCGAL', '-lgmp']
    
    # Create the list of extension modules
    ext_modules = []
    # Generation C++ 
    generation_cpp = [os.path.join('HemeLbSetupTool/Model/Generation', cpp)
                      for cpp in ['Neighbours.cpp',
                                  'BuildCGALPolygon.cpp',
                                  'Block.cpp',
                                  'BlockWriter.cpp',
                                  'BufferPool.cpp',
                                  'GeometryGenerator.cpp',
                                  'GeometryWriter.cpp',
                                  'Domain.cpp',
                                  'Site.cpp',
                                  'InconsistentFluidnessError.cpp',
                                  'Index.cpp',
                                  'CylinderGenerator.cpp',
                                  'PolyDataGenerator.cpp',
                                  'SquareDuctGenerator.cpp',
                                  'Debug.cpp']]
    # HemeLB classes
    hemelb_cpp = [os.path.join(HemeLbDir, cpp)
                  for cpp in ['util/Vector3D.cc',
                              'geometry/SiteData.cc',
                              'lb/lattices/D3Q27.cc',
                              'io/formats/geometry.cc',
                              'io/writers/xdr/XdrFileWriter.cc',
                              'io/writers/xdr/XdrMemWriter.cc',
                              'io/writers/xdr/XdrWriter.cc',
                              'io/writers/Writer.cc']]

    # SWIG wrapper
    swig_cpp = ['HemeLbSetupTool/Model/Generation/Wrap.cpp']
    # Do we need to swig it?
    if not os.path.exists(swig_cpp[0]) or os.path.getmtime(swig_cpp[0]) < os.path.getmtime('HemeLbSetupTool/Model/Generation/Wrap.i'):
        cmd = 'swig -I%s -c++ -python -o HemeLbSetupTool/Model/Generation/Wrap.cpp -outdir HemeLbSetupTool/Model HemeLbSetupTool/Model/Generation/Wrap.i' % HemeLbDir
        print cmd
        os.system(cmd)
    
    generation_ext = Extension('HemeLbSetupTool.Model._Generation',
                               sources=generation_cpp + hemelb_cpp + swig_cpp,
                               extra_compile_args=extra_compile_args,
                               include_dirs=include_dirs,
                               extra_link_args=extra_link_args,
                               library_dirs=library_dirs,
                               libraries=libraries,
                               )
    
    setup(name='HemeLbSetupTool',
          version='1.1',
          author='Rupert Nash',
          author_email='rupert.nash@ucl.ac.uk',
          packages=['HemeLbSetupTool', 'HemeLbSetupTool.Bindings', 'HemeLbSetupTool.Util', 'HemeLbSetupTool.Model', 'HemeLbSetupTool.View', 'HemeLbSetupTool.Controller'],
          scripts=['scripts/hemelb-setup', 'scripts/hemelb-setup-nogui', 'scripts/hemelb-countsites'],
          ext_modules=[generation_ext]
          )
