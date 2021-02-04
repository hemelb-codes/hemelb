# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(parmetis find)
if (PARMETIS_FOUND AND METIS_FOUND)
  message("Parmetis already installed, no need to (download and) build")
  add_custom_target(parmetis)
else()
  message("Parmetis not installed, will build from source")
  find_file(PARMETIS_TARBALL parmetis-4.0.2.tar.gz 
    DOC "Path to download Parmetis (can be url http://)"
    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
    )
  if(NOT PARMETIS_TARBALL)
    message("No parmetis source found, will download.")
    set(PARMETIS_TARBALL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.2.tar.gz 
      CACHE STRING "Path to download Parmetis (can be local file://)" FORCE)
  endif()
  
  set(PARMETIS_CC ${MPI_C_COMPILER} CACHE STRING "MPI Compiler to use for Parmetis, leave blank to let parmetis guess")
  set(PARMETIS_CXX ${MPI_CXX_COMPILER} CACHE STRING "MPI Compiler to use for Parmetis, leave blank to let parmetis guess")
  if(PARMETIS_CC)
    set(PARMETIS_CC_OPTION cc=${PARMETIS_CC})
  endif()
  if(PARMETIS_CXX)
    set(PARMETIS_CXX_OPTION cxx=${PARMETIS_CXX})
  endif()
  ExternalProject_Add(
    parmetis
    INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
    URL ${PARMETIS_TARBALL}
    CONFIGURE_COMMAND make config prefix=<INSTALL_DIR> ${PARMETIS_CC_OPTION} ${PARMETIS_CXX_OPTION} && 
    cd metis && 
    make config prefix=<INSTALL_DIR> ${PARMETIS_CC_OPTION} ${PARMETIS_CXX_OPTION}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS} MAKEFLAGS= #Workaround for bug in parmetis makefile treating switches as targets
    INSTALL_COMMAND make install MAKEFLAGS= && cd metis && make install MAKEFLAGS=
    )
endif()
