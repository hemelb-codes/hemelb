add_library(hemelb_lb
	boundaries/BoundaryComms.cc boundaries/BoundaryValues.cc
	boundaries/iolets/InOutLetCosine.cc boundaries/iolets/InOutLetFile.cc
	boundaries/iolets/InOutLetMultiscale.cc boundaries/iolets/InOutLetVelocityAware.cc
	boundaries/iolets/InOutLetParabolicVelocity.cc
	IncompressibilityChecker.cc
	kernels/momentBasis/DHumieresD3Q15MRTBasis.cc kernels/momentBasis/DHumieresD3Q19MRTBasis.cc
	kernels/rheologyModels/AbstractRheologyModel.cc kernels/rheologyModels/CarreauYasudaRheologyModel.cc 
	kernels/rheologyModels/CassonRheologyModel.cc kernels/rheologyModels/TruncatedPowerLawRheologyModel.cc
	lattices/LatticeInfo.cc lattices/D3Q15.cc lattices/D3Q19.cc lattices/D3Q27.cc lattices/D3Q15i.cc
	MacroscopicPropertyCache.cc SimulationState.cc StabilityTester.cc
	 )
