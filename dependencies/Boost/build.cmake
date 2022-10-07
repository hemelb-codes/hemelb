# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

find_file(BOOST_TARBALL boost_1_75_0.tar.gz
  DOC "Path to download BOOST (can be url http://)"
  PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
  )
if(NOT BOOST_TARBALL)
  message("No boost source found, will download.")
  set(BOOST_TARBALL https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz
    CACHE STRING "Path to download BOOST (can be local file://)" FORCE)
endif()
ExternalProject_Add(
  dep_Boost
  INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PREFIX}
  URL ${BOOST_TARBALL}
  BUILD_COMMAND ""
  INSTALL_COMMAND cp -r <SOURCE_DIR>/boost <INSTALL_DIR>/include
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
)
