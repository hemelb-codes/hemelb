#!/bin/sh
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

rm -rf CMakeCache.txt CMakeFiles setup.cfg Makefile cmake_install.cmake
rm -rf HemeLbSetupTool/Model/Generation_wrap.cpp HemeLbSetupTool/Model/Generation.py HemeLbSetupTool/Model/_Generation.so
rm -rf HemeLbSetupTool/Model/Test_wrap.cpp HemeLbSetupTool/Model/Test.py HemeLbSetupTool/Model/_Test.so
rm -rf build Code
