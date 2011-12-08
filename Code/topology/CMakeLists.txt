if(USE_MULTIMACHINE)
	set(topology MultiMachineNetworkTopology.cc)
else(USE_MULTIMACHINE)
	set(topology SingleMachineNetworkTopology.cc)
endif(USE_MULTIMACHINE)
add_library(hemelb_topology 
	NetworkTopology.cc 
	${topology}
)