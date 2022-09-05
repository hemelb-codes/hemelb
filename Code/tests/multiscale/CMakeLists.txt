# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
add_test_lib(test_multiscale
  MockIntercommunicand.cc
  MockIntercommunicator.cc
  MockIntercommunicatorTests.cc
)

if (HEMELB_BUILD_MULTISCALE)
  add_subdirectory(mpwide)
endif()
