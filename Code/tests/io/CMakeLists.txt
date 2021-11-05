# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/xdr_test_data.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/PathManagerTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/XdrWriterTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/XdrReaderTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/xml.cc
)
