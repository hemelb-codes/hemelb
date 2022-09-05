# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
add_test_lib(test_redblood
  Fixtures.cc

  BendingTests.cc
  BordersTests.cc
  CellArmyTests.cc
  CellIOTests.cc
  CellInserterTests.cc
  CellCellInteractionTests.cc
  CellCellInteractionWithGridTests.cc
  CellForceSpreadTests.cc
  CellTests.cc
  CellVelocityInterpolTests.cc
  DivideConquerTests.cc
  EnergyTests.cc
  FacetTests.cc
  FaderCellTests.cc
  FlowExtensionTests.cc
  GradientKernTests.cc
  GradientTests.cc
  InterpolationTests.cc
  LoadDeformedCellTests.cc
  LoadingTimmMeshTests.cc
  Node2NodeTests.cc
  NodeIntegrationTests.cc
  RedBloodMeshDataIOTests.cc
  RedBloodMeshTests.cc
  RedBloodMeshVTKDataIOTests.cc
  StencilTests.cc
  TopologyTests.cc
  VertexBagTests.cc
  WallCellPairIteratorTests.cc
  WallNodeDnCTests.cc

  buffer/BufferTests.cc
  buffer/ColumnsTests.cc

  parallel/NodeParallelizationTests.cc

  # Longer-running integration tests
  CellIntegrationTests.cc
  FadeInOutIntegrationTests.cc
  SadCellTest.cc
  SadCellIntegrationTests.cc
  MulticellBenchmarkTests.cc
  FedosovValidationTests.cc
)

target_link_libraries(test_redblood PRIVATE VTK::CommonCore)

add_test_lib(test_mpi_redblood
  Fixtures.cc
  parallel/ParallelFixtures.cc
  parallel/ParallelFixtureTests.cc

  parallel/CellParallelizationTests.cc
  parallel/MPISpreadForcesTests.cc
  parallel/MPIIntegrateVelocities.cc
  parallel/MPILockStepTests.cc
  parallel/GraphCommsTests.cc
  parallel/RBCInserterWithPerturbationParallelTests.cc
  parallel/MPIParallelIntegrationTests.cc
)
target_link_libraries(test_mpi_redblood PRIVATE VTK::CommonCore)