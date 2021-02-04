# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(catch2 find)
if (Catch2_FOUND)
  message("Catch2 already installed, no need to download")
  add_custom_target(catch2)
else()
  message("Catch2 not installed, will build from source")

  find_file(CATCH2_TARBALL v2.13.0.tar.gz 
    DOC "Path to download Catch2 (can be url http://)"
    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
    )
  if(NOT CATCH2_TARBALL)
    message("No Catch2 source found, will download.")
    set(CATCH2_TARBALL https://github.com/catchorg/Catch2/archive/v2.13.0.tar.gz
      CACHE STRING "Path to download Catch2 (can be local file://)" FORCE)
  endif()
  ExternalProject_Add(
    catch2
    URL ${CATCH2_TARBALL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${HEMELB_DEPENDENCIES_INSTALL_PATH} -DBUILD_TESTING=OFF
    )
endif()
