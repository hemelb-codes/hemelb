# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.


set(EXTRACTION_SOURCES
  GeometrySelector.cc
  StraightLineGeometrySelector.cc
  IterableDataSource.cc
  PlaneGeometrySelector.cc
  PropertyActor.cc
  PropertyWriter.cc
  WholeGeometrySelector.cc
  LbDataSourceIterator.cc
  GeometrySurfaceSelector.cc
  SurfacePointSelector.cc
  LocalDistributionInput.cc
  
  # reconstructed base class from the original class
  LocalPropertyOutput.cc
  # subclass for original XDR functionality
  XdrPropertyOutput.cc
  )

# new source files for HDF5 support
# only include these if USE_HDF5 is enabled
if(USE_HDF5)
  list(APPEND EXTRACTION_SOURCES
    LocalPropertyHdf5Output.cc
  )
endif()

# use the list of sources to create the hemelb_extraction library
add_library(hemelb_extraction OBJECT
  ${EXTRACTION_SOURCES}
)

# new CMake module for HDF5 support
# only include this if USE_HDF5 is enabled
if(USE_HDF5)
  # include the HDF5 module into hemelb_extraction module's list of dependencies
  # PUBLIC makes the include directories and libraries available to consumers of hemelb_extraction
  target_include_directories(hemelb_extraction PUBLIC
    ${HDF5_INCLUDE_DIRS}
  )

  # link the hemelb_extraction library with HDF5
  target_link_libraries(hemelb_extraction PUBLIC
    HDF5::HDF5
  )
  # define a preprocessor macro to indicate that HDF5 is being used
  target_compile_definitions(hemelb_extraction PUBLIC USE_HDF5)
endif()