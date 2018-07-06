# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import sys
import os.path
import warnings
from setuptools import setup
from setuptools.extension import Extension

def conda_monkey_patch():
    import setuptools.command.easy_install
    BaseScriptWriter = setuptools.command.easy_install.ScriptWriter
    
    class CondaOsxScriptWriter(BaseScriptWriter):
        @classmethod
        def _get_script_args(cls, type_, name, header, script_text):
            
            if type_ == 'gui':
                warnings.warn('Using custom hack to work around Conda on macOS for scripts that launch a GUI', UserWarning)
                old_exe = cls.command_spec_class.best().from_param(None)[0]
                dirname, python = os.path.split(old_exe)
                new_exe = os.path.join(dirname, 'python.app')
                assert os.path.exists(new_exe)
                header = header.replace(old_exe, new_exe)

            yield (name, header + script_text)
        pass
    
    setuptools.command.easy_install.BaseScriptWriter = BaseScriptWriter
    setuptools.command.easy_install.ScriptWriter = CondaOsxScriptWriter
    
    return

def on_conda():
    '''Try to figure out if we running on a conda interpreter by looking
    for the conda command line tool in the same directory as this
    interpreter.
    '''
    bindir, python = os.path.split(sys.executable)
    return os.path.exists(os.path.join(bindir, 'conda'))
    
    
if sys.platform == 'darwin':
    # Python thinks it's so smart and sets the
    # MACOSX_DEPLOYMENT_TARGET environment variable that messes around
    # with what features of the compiler and C++ std lib are
    # available. Set this to use your current one.
    import platform
    release, versioninfo, machine = platform.mac_ver()
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = release

    if on_conda():
        conda_monkey_patch()
    

lib_src = '''
../../Code/util/Vector3D.cc
HemeLbSetupTool/Model/Generation/FloodFill.cpp
HemeLbSetupTool/Model/Generation/H5.cpp
HemeLbSetupTool/Model/Generation/MkCgalMesh.cpp
HemeLbSetupTool/Model/Generation/Neighbours.cpp
HemeLbSetupTool/Model/Generation/PolyDataGenerator.cpp
HemeLbSetupTool/Model/Generation/SectionTree.cpp
HemeLbSetupTool/Model/Generation/SectionTreeBuilder.cpp
HemeLbSetupTool/Model/Generation/SegmentFactory.cpp
HemeLbSetupTool/Model/Generation/SurfaceVoxeliser.cpp
HemeLbSetupTool/Model/Generation/TriTree.cpp
HemeLbSetupTool/Model/Generation/TriangleSorter.cpp
HemeLbSetupTool/Model/Generation/Vector.cpp
'''.strip().split('\n')


generation_src = lib_src + ['HemeLbSetupTool/Model/Generation.i']
test_src = lib_src + ['HemeLbSetupTool/Model/Generation/test.cpp',
                          'HemeLbSetupTool/Model/Test.i']


main_libs = ['boost_system', 'boost_filesystem', 'hdf5', 'CGAL', 'gmp', 'mpfr']

generation_ext = Extension(
    'HemeLbSetupTool.Model._Generation',
    sources=generation_src,
    libraries=main_libs,
    extra_compile_args=['--std=c++11']
    )
test_ext = Extension(
    'HemeLbSetupTool.Model._Test',
    sources=test_src,
    libraries=main_libs+['cppunit'],
    extra_compile_args=['--std=c++11']
    )

setup(
    name='HemeLbSetupTool',
    version='2.0',
    author='Rupert Nash',
    author_email='r.nash@epcc.ed.ac.uk',
    packages=[
        'HemeLbSetupTool',
        'HemeLbSetupTool.Bindings',
        'HemeLbSetupTool.Util',
        'HemeLbSetupTool.Model',
        'HemeLbSetupTool.View',
        'HemeLbSetupTool.Controller',
        'HemeLbSetupTool.scripts'
        ],
        # Define entry points instead of scripts
        # https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
    entry_points={
        'console_scripts': [
            'hemelb-setup-nogui=HemeLbSetupTool.scripts.setup_cli:main',
            'hemelb-countsites=HemeLbSetupTool.scripts.countsites:main',
            'upgrade-profile=HemeLbSetupTool.scripts.upgrade_profile:main'
            ],
        'gui_scripts': [
            'hemelb-setup=HemeLbSetupTool.scripts.setup_gui:main'
            ]
        },
    ext_modules=[generation_ext, test_ext]
    )
