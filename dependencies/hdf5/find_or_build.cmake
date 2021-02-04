# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(hdf5 find)

if (HDF5_FOUND AND HDF5_IS_PARALLEL)
  message("HDF5 with MPI support found, no need to download")
else()
  message("HDF5 with MPI support not installed, will build from source")

  find_file(HDF5_TARBALL hdf5-1.8.15-patch1.tar.gz
    DOC "Path to download HDF5 (can be url http://)"
    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
    )
  if (NOT HDF5_TARBALL)
    message("No HDF5 source found, will download")
    set(HDF5_TARBALL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.15-patch1/src/hdf5-1.8.15-patch1.tar.gz
      CACHE STRING "Path to download HDF5 (can be local file://)" FORCE)
  endif()
  
  ExternalProject_Add(
    hdf5
    INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
    URL ${HDF5_TARBALL}
    CONFIGURE_COMMAND cmake <SOURCE_DIR>
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
    -DBUILD_SHARED_LIBS=ON
    -DBUILD_TESTING=OFF
    -DHDF5_BUILD_CPP_LIB=OFF
    -DHDF5_BUILD_EXAMPLES=OFF
    -DHDF5_ENABLE_PARALLEL=ON
    -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
    BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
    )
endif()
