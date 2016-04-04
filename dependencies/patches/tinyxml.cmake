
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
cmake_minimum_required(VERSION 2.6)

project(TinyXml)
OPTION(TIXML_USE_STL "Use STL with TIXML" ON)
if(TIXML_USE_STL)
	add_definitions(-DTIXML_USE_STL)
endif(TIXML_USE_STL)
add_library(
	tinyxml
	tinyxml.cpp
	tinystr.cpp
	tinyxmlerror.cpp
	tinyxmlparser.cpp
)

INSTALL( FILES tinyxml.h tinystr.h DESTINATION include )
INSTALL( TARGETS tinyxml ARCHIVE DESTINATION lib )
