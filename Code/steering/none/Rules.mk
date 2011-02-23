include $(MK)/header.mk

SRCS := SteeringComponent.cc \
        ImageSendComponent.cc \
        ClientConnection.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
