# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_geometry OBJECT
  GmyReadResult.cc
  BlockTraverser.cc
  GeometryReader.cc needs/Needs.cc
  LookupTree.cc
  Domain.cc FieldData.cc
  SiteDataBare.cc SiteData.cc
  SiteTraverser.cc VolumeTraverser.cc Block.cc
  decomposition/BasicDecomposition.cc
  decomposition/OptimisedDecomposition.cc
        neighbouring/NeighbouringDomain.cc
  neighbouring/NeighbouringDataManager.cc
  neighbouring/RequiredSiteInformation.cc
  )
target_link_libraries(hemelb_geometry PRIVATE ParMETIS::ParMETIS ZLIB::ZLIB)
target_link_libraries(hemelb_geometry PUBLIC Boost::boost)
configure_file (
  "${PROJECT_SOURCE_DIR}/geometry/decomposition/DecompositionWeights.h.in"
  "${PROJECT_BINARY_DIR}/geometry/decomposition/DecompositionWeights.h"
  )

configure_file (
  "${PROJECT_SOURCE_DIR}/geometry/ParmetisForward.h.in"
  "${PROJECT_BINARY_DIR}/geometry/ParmetisForward.h"
  )
