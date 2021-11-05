# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(tinyxml find)
if (TINYXML_FOUND)
  message("TinyXML already installed, no need to download")
  add_custom_target(tinyxml)
else()
  message("TinyXML not installed, will build from source")
  find_file(TINYXML_TARBALL tinyxml_2_6_2.tar.gz 
    DOC "Path to download TinyXML (can be url http://)"
    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
    )
  if(NOT TINYXML_TARBALL)
    message("No tinyxml source found, will download.")
    set(TINYXML_TARBALL http://sourceforge.net/projects/tinyxml/files/tinyxml/2.6.2/tinyxml_2_6_2.tar.gz
      CACHE STRING "Path to download TinyXML (can be local file://)" FORCE)
  endif()
  ExternalProject_Add(
    tinyxml
    INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
    URL ${TINYXML_TARBALL}
    CONFIGURE_COMMAND cmake <SOURCE_DIR>
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DTIXML_USE_STL=ON
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
    PATCH_COMMAND cp ${HEMELB_DEPENDENCIES_PATH}/patches/tinyxml.cmake CMakeLists.txt
    BUILD_COMMAND make
    )
endif()
