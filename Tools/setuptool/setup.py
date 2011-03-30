import os.path
import sys
import platform
import re

import numpy
import vtk

from distutils.core import setup
# from distutils.extension import Extension

from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension

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
        # Get the second half, showing where the symbol found
        name, pathAndStuff = line.split('=>')
        
        # Split on whitespace
        words = pathAndStuff.split()
        try:
            # First word is the full path of the library
            libPath = words[0]
        except IndexError:
            continue
        try:
            dir, base = os.path.split(libPath)
            # Do a string split rather than os.path.splitext as the latter splits on the last dot
            base, rest = base.split('.')[0]
            if base == 'libvtkCommon':
                return dir
            
        except ValueError:
            pass
        continue
    # Searched them all without finding a match!
    raise ValueError("Can't figure out the VTK include dir")

def LibToInclude(vtkLibDir):
    libDir, vtk = os.path.split(vtkLibDir)
    prefix, lib = os.path.split(libDir)
    includeDir = os.path.join(prefix, 'include')
    vtkIncludeDir = os.path.join(includeDir, vtk)
    return vtkIncludeDir

def GetVtkLibDir():
    aVtkSharedLibrary = vtk.libvtkCommonPython.__file__
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

if __name__ == "__main__":
    # numpy, vtk
    vtkLibDir = GetVtkLibDir()
    
    include_dirs = [numpy.get_numpy_include(), LibToInclude(vtkLibDir)]
    libraries = []
    library_dirs = []
    extra_compile_args = GetVtkCompileFlags(vtkLibDir)
    extra_link_args = []

    # Create the list of extension modules
    ext_modules = []
    for pyx in ['HemeLbSetupTool/Model/HitList.pyx',
                'HemeLbSetupTool/Model/OutputGenerationHelpers.pyx',
                'HemeLbSetupTool/Model/vtkHelp.pyx']:
        # Work out the module name
        base, ignored = os.path.splitext(pyx)
        modname = base.replace('/', '.')

        ext = Extension(modname,
                        [pyx],
                        extra_compile_args=extra_compile_args,
                        include_dirs=include_dirs,
                        extra_link_args=extra_link_args,
                        library_dirs=library_dirs,
                        libraries=libraries,
                        pyrex_cplus=True,
                        )
        ext_modules.append(ext)


    setup(name='HemeLbSetupTool',
          version='1.0',
          author='Rupert Nash',
          author_email='rupert.nash@ucl.ac.uk',
          packages=['HemeLbSetupTool', 'HemeLbSetupTool.Bindings', 'HemeLbSetupTool.Util', 'HemeLbSetupTool.Model', 'HemeLbSetupTool.View', 'HemeLbSetupTool.Controller'],
          scripts=['scripts/setuptool'],
          # Swap in the Cython builder
          cmdclass={'build_ext': build_ext},
          ext_modules=ext_modules
          )
