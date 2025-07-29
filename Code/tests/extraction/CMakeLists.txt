# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
add_test_lib(test_extraction
  GeometrySelectorTests.cc
  LocalPropertyOutputTests.cc
  LocalPropertyHdf5OutputTests.cc
  )

  # If HDF5 support is enabled in the main build,
  # also enable it for compiling these tests.
  if(USE_HDF5)
  target_compile_definitions(test_extraction
    PUBLIC USE_HDF5
  )
  endif()