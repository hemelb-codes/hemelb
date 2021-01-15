# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

if (HEMELB_BUILD_RBC)
  # Build full RBC code
  hemelb_dependency(hdf5 use)
  hemelb_dependency(vtk use)

  add_library(hemelb_redblood
    RBCConfig.cc
    Mesh.cc CellBase.cc Cell.cc Interpolation.cc
    CellCell.cc FlowExtension.cc FaderCell.cc io.cc RBCInserter.cc
    VertexBag.cc Borders.cc
    WallCellPairIterator.cc
    VTKError.cc
    buffer/Buffer.cc buffer/Columns.cc
    parallel/NodeCharacterizer.cc parallel/CellParallelization.cc parallel/SpreadForces.cc
    parallel/IntegrateVelocities.cc parallel/GraphBasedCommunication.cc
    )

  hemelb_add_target_dependency_hdf5(hemelb_redblood)
  hemelb_add_target_dependency_vtk(hemelb_redblood)
  hemelb_add_target_dependency_boost(hemelb_redblood)
else()
  # Build minimal interface
  add_library(hemelb_redblood
    RBCConfig.cc
    )
endif()
