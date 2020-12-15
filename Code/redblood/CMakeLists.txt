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
    Mesh.cc MeshIO.cc
    CellBase.cc Cell.cc CellEnergy.cc Facet.cc
    Interpolation.cc
    CellCell.cc FlowExtension.cc FaderCell.cc xmlIO.cc RBCInserter.cc
    VertexBag.cc Borders.cc
    WallCellPairIterator.cc
    VTKError.cc
    buffer/Buffer.cc buffer/Columns.cc
    parallel/NodeCharacterizer.cc parallel/CellParallelization.cc parallel/SpreadForces.cc
    parallel/IntegrateVelocities.cc parallel/GraphBasedCommunication.cc
    )

  target_link_libraries(hemelb_redblood PRIVATE
    hdf5::hdf5_hl
    VTK::CommonCore VTK::CommonDataModel VTK::IOXML VTK::FiltersCore
    Boost::boost)
else()
  # Build minimal interface
  add_library(hemelb_redblood
    RBCConfig.cc
    )
endif()
