#------Parmetis  ------------
find_package(Parmetis REQUIRED)
include_directories(SYSTEM ${PARMETIS_INCLUDE_DIRS})

#------TinyXML ----------------
find_package(TinyXML REQUIRED)
OPTION(TIXML_USE_STL "Use STL with TIXML" ON)
if(TIXML_USE_STL)
	add_definitions(-DTIXML_USE_STL)
endif()
include_directories(SYSTEM ${TINYXML_INCLUDE_DIR})
#------BOOST ------------------
find_package(Boost 1.54 REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
add_definitions(-DHEMELB_USE_BOOST)
#------CTemplate ----------------
find_package(CTemplate REQUIRED)
include_directories(SYSTEM ${CTEMPLATE_INCLUDE_DIR})
if(HEMELB_BUILD_MULTISCALE)
  #------MPWide ----------------
  find_package(MPWide REQUIRED)
  include_directories(SYSTEM ${MPWide_INCLUDE_DIR})
  add_definitions(-DHEMELB_BUILD_MULTISCALE)
endif()

#------zlib ----------------
find_package(ZLIB REQUIRED)
include_directories(SYSTEM ${ZLIB_INCLUDE_DIR})

#------VTK ----------------
find_package(VTK REQUIRED)
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
