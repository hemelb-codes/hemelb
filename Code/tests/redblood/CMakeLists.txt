# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/Fixtures.cc

  ${CMAKE_CURRENT_SOURCE_DIR}/BendingTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/CellCellInteractionTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/CellCellInteractionWithGridTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/CellTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/EnergyTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/FacetTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/GradientKernTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/GradientTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RedBloodMeshDataIOTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RedBloodMeshTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/StencilTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/TopologyTests.cc
  )
