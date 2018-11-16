#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os
import re
import sys
import shutil
from contextlib import contextmanager
import re

py_coding_re = re.compile('#.*coding[=:]\s*([-\w.]+)')

commentchars={
    '.applescript': '--',
    '.cc' : '//',
    '.cmake': '#',
    '.cpp' : '//',
    '.sh':'#',
    '.i' : '//',
    '.hpp' : '//',
    '.h' : '//',
    '.py' : '#',
    '.pyx' : '#',
    '.pxd' : '#',
    'CMakeLists.txt' : '#',
    '.yml' : '#',
    'README': '',
    }

oldcopyblob = """
Copyright (C) University College London, 2007-2012, all rights reserved.

This file is part of HemeLB and is CONFIDENTIAL. You may not work
with, install, use, duplicate, modify, redistribute or share this
file, or any part thereof, other than as allowed by any agreement
specifically made by you with University College London.
"""

newcopyblob = """This file is part of HemeLB and is Copyright (C)
the HemeLB team and/or their institutions, as detailed in the
file AUTHORS. This software is provided under the terms of the
license in the file LICENSE."""

old_commented_cr = {}
new_commented_cr = {}
for ext in commentchars:
    old_commented_cr[ext] = "\n".join([commentchars[ext] + ' ' + line for line in oldcopyblob.split('\n')])
    new_commented_cr[ext] = "\n".join([commentchars[ext] + ' ' + line for line in newcopyblob.split('\n')])
    continue

@contextmanager
def WorkingDirectory(dst_dir):
    current_dir = os.getcwd()
    os.chdir(dst_dir)
    try:
        yield dst_dir
    finally:
        os.chdir(current_dir)


def Walk(codeDir):
    for dirpath, dirnames, filenames in os.walk(codeDir):
        if 'build' in dirnames:
            dirnames.remove('build')
            pass

        for name in filenames:
            base, ext = os.path.splitext(name)
            if ext == '.in':
                ext = os.path.splitext(base)[1]
            if ext in commentchars:
                yield os.path.normpath(os.path.join(dirpath, name)), ext
            elif name == 'CMakeLists.txt':
                yield os.path.normpath(os.path.join(dirpath, name)), name

def SplitAll(path):
    head, tail = os.path.split(path)
    if head == '':
        ans = []
    else:
        ans = SplitAll(head)

    ans.append(tail)
    return ans

old_cr_lines = map(str.strip, oldcopyblob.split('\n'))
nlines = len(old_cr_lines)
def StripOldCopyright(text, comment):
    header = text.split('\n', nlines)
    main_text = header.pop()

    assert len(header) == nlines
    for cr in old_cr_lines:
        line = header.pop(0)
        assert line.startswith(comment)
        assert line[len(comment):].strip() == cr
    return main_text

not_ours = set("""Tools/setuptool/HemeLbSetupTool/Model/Generation/Point_inside_polyhedron_3.h
Tools/setuptool/HemeLbSetupTool/Model/Generation/Triangle_3_Ray_3_do_intersect.h""".split('\n'))

no_cr = set("""Code/debug/linux/launchGdbs.sh
Code/debug/linux/launchGnomeTerminal.sh
Code/debug/linux/launchKonsole.sh
Code/extraction/ArbitrarySiteListIterableDataSource.cc
Code/extraction/SurfacePointSelector.cc
Code/extraction/SurfacePointSelector.h
Code/geometry/ParmetisHeader.h
Code/lb/iolets/InOutLet.cc
Code/lb/iolets/InOutLetFileVelocity.cc
Code/lb/iolets/InOutLetFileVelocity.h
Code/lb/iolets/InOutLetParabolicVelocity.cc
Code/lb/iolets/InOutLetVelocity.cc
Code/lb/iolets/InOutLetWomersleyVelocity.cc
Code/lb/lattices/Lattices.h
Code/lb/streamers/NashZerothOrderPressureDelegate.h
Code/lb/streamers/NashZerothOrderPressureIolet.h
Code/unittests/util/UnitConverterTests.h
Code/util/FlatMap.h
deploy/myfabfile.py
Scripts/ParameterChooser.py
Tools/analysis/test/__init__.py
Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py
Tools/hemeTools/simconf/__init__.py
Tools/hemeTools/simconf/HemeLbParameters.py
Tools/hemeTools/simconf/OneInletHelper.py
Tools/hemeTools/simconf/PoiseuilleResistance.py
Tools/hemeTools/simconf/simplify.py
Tools/hemeTools/simconf/VesselNetwork.py
Tools/hemeTools/utils/__init__.py
Tools/hemeTools/utils/default_property.py
Tools/hemeTools/utils/IterPairs.py
Tools/hemeTools/utils/memo.py
Tools/hemeTools/utils/SortedVector.py
Tools/hemeTools/utils/utils.py
Tools/hemeTools/utils/test/__init__.py
Tools/hemeTools/utils/test/test_utils.py
Tools/hemeTools/vtkhelp/__init__.py
Tools/hemeTools/vtkhelp/cells.py
Tools/setuptool/HemeLbSetupTool/Model/Generation/CGALtypedef.h
Tools/setuptool/InletProcessing/__init__.py
Tools/setuptool/InletProcessing/PlotWeightsFile.py
Tools/setuptool/test/TestHemeLBSetupTool/Model/test_OutputGeneration.py
Tools/setuptool/test/TestHemeLBSetupTool/Model/fixtures/__init__.py
Tools/setuptool/test/TestHemeLBSetupTool/Model/fixtures/__main__.py
Tools/setuptool/test/TestHemeLBSetupTool/Model/fixtures/generate.py
Tools/setuptool/test/TestHemeLBSetupTool/Model/fixtures/load.py
Tools/setuptool/test/TestHemeLBSetupTool/Util/test_Observer.py
CMakeLists.txt
Code/CMakeLists.txt
Code/cmake/build_environment.cmake
Code/cmake/dependencies.cmake
Code/cmake/gnu_bug.cmake
Code/cmake/mountain_lion_scandir.cmake
Code/cmake/mpi.cmake
Code/cmake/platform_checks.cmake
Code/colloids/CMakeLists.txt
Code/configuration/CMakeLists.txt
Code/debug/CMakeLists.txt
Code/debug/OSX/MPIdebug.applescript
Code/extraction/CMakeLists.txt
Code/functionaltests/cpptests/CMakeLists.txt
Code/geometry/CMakeLists.txt
Code/io/CMakeLists.txt
Code/lb/CMakeLists.txt
Code/log/CMakeLists.txt
Code/multiscale/CMakeLists.txt
Code/net/CMakeLists.txt
Code/reporting/CMakeLists.txt
Code/resources/CMakeLists.txt
Code/steering/CMakeLists.txt
Code/unittests/CMakeLists.txt
Code/util/CMakeLists.txt
Code/vis/CMakeLists.txt
dependencies/CMakeLists.txt
Tools/hemeTools/utils/cxdr.pxd
Tools/hemeTools/utils/xdr.pxd
Tools/hemeTools/utils/xdr.pyx
dependencies/Modules/FindCPPUnit.cmake
dependencies/Modules/FindCTemplate.cmake
dependencies/Modules/FindMetis.cmake
dependencies/Modules/FindMPI.cmake
dependencies/Modules/FindMPWide.cmake
dependencies/Modules/FindParmetis.cmake
dependencies/Modules/FindTinyXML.cmake
dependencies/patches/tinyxml.cmake""".split('\n'))

if __name__ == '__main__':
    script = os.path.abspath(sys.argv[0])
    scriptDir = os.path.dirname(script)
    codeDir = os.path.normpath(os.path.join(scriptDir, os.pardir))

    with WorkingDirectory(codeDir):
        for sourceFile, ext in Walk('.'):
            if sourceFile in not_ours:
                continue

            print sourceFile

            oldtext = open(sourceFile).read()

            if oldtext[0:2] == '#!':
                envline, oldtext = oldtext.split('\n', 1)
            else:
                envline = ''
                pass

            if ext == '.py':
                # Python files may have a file encoding specified in the first two lines
                if py_coding_re.match(oldtext):
                    codingline, oldtext = oldtext.split('\n', 1)
                    if envline:
                        envline += '\n' + codingline
                    else:
                        envlince = codingline
                        pass
                    pass
                pass

            new_cr = new_commented_cr[ext]
            if sourceFile in no_cr:
                main = oldtext
            else:
                try:
                    main = StripOldCopyright(oldtext, commentchars[ext])
                except AssertionError:
                    sys.stderr.write(sourceFile + '\n')
                    continue
                pass

            with open(sourceFile, 'w') as replacement:
                replacement.write(envline + '\n')
                replacement.write(new_cr + '\n')
                replacement.write(main)
