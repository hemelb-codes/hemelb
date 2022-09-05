# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(test-helpers OBJECT
  FolderTestFixture.cc
  FourCubeBasedTestFixture.cc
  FourCubeLatticeData.cc
  HasCommsTestFixture.cc
  LabelledRequest.cc
  LatticeDataAccess.cc
  MockNetHelper.cc
  OneInOneOutSimConfig.cc
  RandomSourceTests.cc
  RecordingNet.cc
)

target_link_libraries(test-helpers
  PRIVATE
  Catch2::Catch2
  TinyXML::TinyXML
)
