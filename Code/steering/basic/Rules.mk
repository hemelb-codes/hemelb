include $(MK)/header.mk

SRCS := HttpPost.cc \
        Network.cc \
        SimulationParameters.cc \
        Threadable.cc \
        ImageSendComponent.cc \
        ClientConnection.cc \
        SteeringComponent.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
