# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

file(GLOB_RECURSE HEMELB_CU_FILES ${CMAKE_SOURCE_DIR}/*.cc)
set_source_files_properties( ${HEMELB_CU_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ )

cuda_add_library(hemelb_lb
  iolets/BoundaryCommunicator.cc iolets/BoundaryComms.cc iolets/BoundaryValues.cc
  iolets/InOutLet.cc
  iolets/InOutLetCosine.cc iolets/InOutLetFile.cc
  iolets/InOutLetMultiscale.cc
  iolets/InOutLetVelocity.cc
  iolets/InOutLetParabolicVelocity.cc iolets/InOutLetWomersleyVelocity.cc iolets/InOutLetFileVelocity.cc
  IncompressibilityChecker.cc
  kernels/momentBasis/DHumieresD3Q15MRTBasis.cc kernels/momentBasis/DHumieresD3Q19MRTBasis.cc
  kernels/rheologyModels/AbstractRheologyModel.cc kernels/rheologyModels/CarreauYasudaRheologyModel.cc 
  kernels/rheologyModels/CassonRheologyModel.cc kernels/rheologyModels/TruncatedPowerLawRheologyModel.cc
  lattices/LatticeInfo.cc lattices/D3Q15.cc lattices/D3Q19.cc lattices/D3Q27.cc lattices/D3Q15i.cc
  MacroscopicPropertyCache.cc SimulationState.cc StabilityTester.cc
  streamers/StreamerTypeFactory.cu
  )
