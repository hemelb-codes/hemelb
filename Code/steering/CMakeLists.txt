if(HEMELB_STEERING_LIB MATCHES [Bb]asic)
	set(steerers 
		basic/ClientConnection.cc
		basic/HttpPost.cc
		basic/ImageSendComponent.cc
		basic/Network.cc
		basic/SimulationParameters.cc
		basic/SteeringComponent.cc
	)
	add_definitions(-DHEMELB_STEERING_LIB=basic)
elseif(HEMELB_STEERING_LIB MATCHES [Nn]one)
	set(steerers
		none/ClientConnection.cc
		none/ImageSendComponent.cc
		none/Network.cc
		none/SteeringComponent.cc
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
		basic/SteeringComponent.cc
	)
	add_definitions(-DHEMELB_STEERING_LIB=basic)
endif()

add_library(hemelb_steering
	#common/Steerer.cc # Not used in old nrmake build either -- TODO find out why
	common/SteeringComponent.cc
	#common/Tags.cc
	${steerers}
)