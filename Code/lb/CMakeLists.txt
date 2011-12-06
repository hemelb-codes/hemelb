add_library(hemelb_lb
	boundaries/BoundaryComms.cc boundaries/BoundaryValues.cc
	boundaries/iolets/InOutLet.cc boundaries/iolets/InOutLetCosine.cc boundaries/iolets/InOutLetFile.cc
	EntropyTester.cc HFunction.cc IncompressibilityChecker.cc io.cc
	kernels/rheologyModels/AbstractRheologyModel.cc kernels/rheologyModels/CarreauYasudaRheologyModel.cc 
	kernels/rheologyModels/CassonRheologyModel.cc kernels/rheologyModels/TruncatedPowerLawRheologyModel.cc
	lb.cc SimulationState.cc StabilityTester.cc
	 )