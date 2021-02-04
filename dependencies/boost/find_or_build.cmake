# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(boost find)
if (Boost_FOUND)
  message("Boost >= 1.54 already installed, no need to download")
  add_custom_target(boost)
else()
  message("Boost not installed, will build from source")
  find_file(BOOST_TARBALL boost_1_60_0.tar.gz
    DOC "Path to download BOOST (can be url http://)"
    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
    )
  if(NOT BOOST_TARBALL)
    message("No boost source found, will download.")
    set(BOOST_TARBALL http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
      CACHE STRING "Path to download BOOST (can be local file://)" FORCE)
  endif()
  ExternalProject_Add(
    boost
    INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
    URL ${BOOST_TARBALL}
    BUILD_COMMAND ""
    INSTALL_COMMAND cp -r <SOURCE_DIR>/boost <INSTALL_DIR>/include
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    )
endif()
