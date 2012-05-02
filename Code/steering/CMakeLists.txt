set(HEMELB_STEERING_HOST "CCS" CACHE STRING "Use a default host suffix for steering? (CCS, NGS2Leeds, NGS2Manchester, LONI, NCSA or blank)")
if(HEMELB_STEERING_HOST)
	add_definitions("-D${HEMELB_STEERING_HOST}")
endif()

if(HEMELB_STEERING_LIB MATCHES [Bb]asic)
	set(steerers 
		basic/ClientConnection.cc
		basic/HttpPost.cc
		basic/ImageSendComponent.cc
		basic/Network.cc
		basic/SimulationParameters.cc
		basic/SteeringComponentB.cc
	)
	add_definitions(-DHEMELB_STEERING_LIB=basic)
elseif(HEMELB_STEERING_LIB MATCHES [Nn]one)
	set(steerers
		none/ClientConnection.cc
		none/ImageSendComponent.cc
		none/Network.cc
		none/SteeringComponentN.cc
	)
	add_definitions(-DHEMELB_STEERING_LIB=none)
else()
	message("Unrecognised steering mode, using basic")
	set(steerers 
		basic/ClientConnection.cc
		basic/HttpPost.cc
		basic/ImageSendComponent.cc
		basic/Network.cc
		basic/SimulationParameters.cc
		basic/SteeringComponentB.cc
	)
	add_definitions(-DHEMELB_STEERING_LIB=basic)
endif()

add_library(hemelb_steering
	#common/Steerer.cc # Not used in old nrmake build either -- TODO find out why
	common/SteeringComponentC.cc
	#common/Tags.cc
	${steerers}
)
