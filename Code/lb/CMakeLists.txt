
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
add_library(hemelb_lb
	iolets/BoundaryCommunicator.cc iolets/BoundaryComms.cc iolets/BoundaryValues.cc stents/BoundaryValues.cc
	iolets/InOutLet.cc
	iolets/InOutLetCosine.cc iolets/InOutLetFile.cc
	iolets/InOutLetMultiscale.cc
    iolets/InOutLetVelocity.cc
	iolets/InOutLetParabolicVelocity.cc iolets/InOutLetWomersleyVelocity.cc iolets/InOutLetFileVelocity.cc
        stents/Stent.cc
        stents/StentConstant.cc stents/StentFile.cc
        stents/StentFlux.cc
        stents/StentConstantFlux.cc stents/StentFileFlux.cc
	IncompressibilityChecker.cc
	kernels/momentBasis/DHumieresD3Q15MRTBasis.cc kernels/momentBasis/DHumieresD3Q19MRTBasis.cc
	kernels/rheologyModels/AbstractRheologyModel.cc kernels/rheologyModels/CarreauYasudaRheologyModel.cc 
	kernels/rheologyModels/CassonRheologyModel.cc kernels/rheologyModels/TruncatedPowerLawRheologyModel.cc
	lattices/LatticeInfo.cc lattices/D3Q15.cc lattices/D3Q19.cc lattices/D3Q27.cc lattices/D3Q15i.cc
        lattices/D3Q15ad.cc lattices/D3Q19ad.cc lattices/D3Q27ad.cc lattices/D3Q7ad.cc
	MacroscopicPropertyCache.cc SimulationState.cc StabilityTester.cc
	 )
