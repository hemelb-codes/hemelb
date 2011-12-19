if(HEMELB_USE_MULTIMACHINE)
	set(topology MultiMachineNetworkTopology.cc)
else(HEMELB_USE_MULTIMACHINE)
	set(topology SingleMachineNetworkTopology.cc)
endif(HEMELB_USE_MULTIMACHINE)
add_library(hemelb_topology 
	NetworkTopology.cc 
	${topology}
)
