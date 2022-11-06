# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

find_file(CATCH2_TARBALL v2.13.9.tar.gz
  DOC "Path to download Catch2 (can be url http://)"
  PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
  )
if(NOT CATCH2_TARBALL)
  message("No Catch2 source found, will download.")
  set(CATCH2_TARBALL https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.9.tar.gz
    CACHE STRING "Path to download Catch2 (can be local file://)" FORCE)
endif()
ExternalProject_Add(
  dep_Catch2
  URL ${CATCH2_TARBALL}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${HEMELB_DEPENDENCIES_INSTALL_PREFIX} -DBUILD_TESTING=OFF
  )
