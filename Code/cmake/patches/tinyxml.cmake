cmake_minimum_required(VERSION 2.6)

project(TinyXml)
add_definitions(-DTIXML_USE_STL)
add_library(
	tinyxml
	tinyxml.cpp
	tinystr.cpp
	tinyxmlerror.cpp
	tinyxmlparser.cpp
)

INSTALL( FILES tinyxml.h tinystr.h DESTINATION include )
INSTALL( TARGETS tinyxml ARCHIVE DESTINATION lib )
