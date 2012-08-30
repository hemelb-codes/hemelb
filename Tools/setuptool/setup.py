import os.path
import sys
import platform
import re

import numpy
import vtk
try:
    import wxversion
    wxversion.select('2.8')
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

def DarwinGrep(results):
    lines = results.split('\n')
    # Ditch first line
    lines.pop(0)
    
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

def LinuxGrep(results):
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
        aVtkSharedLibrary=vtk.vtkCommonPython.__file__
    osName = platform.system()
    if osName == 'Darwin':
        sharedLibCmd = 'otool -L %s'
        grep = DarwinGrep
    elif osName == 'Linux':
        sharedLibCmd = 'ldd %s'
        grep = LinuxGrep
    else:
        raise ValueError("Don't know how to determing VTK path on OS '%s'" % osName)

    results = os.popen(sharedLibCmd % aVtkSharedLibrary).read()
    return grep(results)

def GetBoostDir(hemeLbDir):
    boostDir = os.path.join(hemeLbDir, '../dependencies/include/') 
    return boostDir

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

def GetHemeLbCompileFlags():
    osName = platform.system()
    if osName == 'Darwin':
        return ['-DHEMELB_CFG_ON_BSD']
    
    return []
    
if __name__ == "__main__":
    # numpy, vtk
    vtkLibDir = GetVtkLibDir()
    HemeLbDir = os.path.abspath('../../Code')
    BoostDir = GetBoostDir(HemeLbDir)
    include_dirs = [ LibToInclude(vtkLibDir), HemeLbDir, BoostDir]
    libraries = []
    library_dirs = []
    extra_compile_args = GetVtkCompileFlags(vtkLibDir) + GetHemeLbCompileFlags()
    extra_link_args = []

    # Create the list of extension modules
    ext_modules = []
    # Generation C++ 
    generation_cpp = [os.path.join('HemeLbSetupTool/Model/Generation', cpp)
                      for cpp in ['Neighbours.cpp',
                                  'Block.cpp',
                                  'BlockWriter.cpp',
                                  'GeometryGenerator.cpp',
                                  'GeometryWriter.cpp',
                                  'Domain.cpp',
                                  'Site.cpp',
                                  'InconsistentFluidnessError.cpp',
                                  'Index.cpp',
                                  'CylinderGenerator.cpp',
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
