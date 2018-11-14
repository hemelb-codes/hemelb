# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_geometry STATIC
  BlockTraverser.cc
  BlockTraverserWithVisitedBlockTracker.cc
  GeometryReader.cc
  needs/Needs.cc
  LatticeData.cc
  LatticeData.cu
  SiteDataBare.cu
  SiteData.cc
  SiteTraverser.cc
  VolumeTraverser.cc
  Block.cc
  decomposition/BasicDecomposition.cc
  decomposition/OptimisedDecomposition.cc
  neighbouring/NeighbouringLatticeData.cc
  neighbouring/NeighbouringDataManager.cc
  neighbouring/RequiredSiteInformation.cc
)
set_property(TARGET hemelb_geometry PROPERTY CUDA_SEPARABLE_COMPILATION ON)

hemelb_add_target_dependency_parmetis(hemelb_geometry)
hemelb_add_target_dependency_zlib(hemelb_geometry)

configure_file (
  "${PROJECT_SOURCE_DIR}/geometry/decomposition/DecompositionWeights.h.in"
  "${PROJECT_BINARY_DIR}/geometry/decomposition/DecompositionWeights.h"
)

configure_file (
  "${PROJECT_SOURCE_DIR}/geometry/ParmetisForward.h.in"
  "${PROJECT_BINARY_DIR}/geometry/ParmetisForward.h"
)
