# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
add_test_lib(test_lb
  BroadcastMocks.cc
  CollisionTests.cc
  IncompressibilityCheckerTests.cc
  KernelTests.cc
  LatticeTests.cc
  RheologyModelTests.cc
  StreamerTests.cc
  VirtualSiteIoletStreamerTests.cc
  GuoForcingTests.cc
  )
add_subdirectory(iolets)
target_link_libraries(test_lb PRIVATE test_iolets)
