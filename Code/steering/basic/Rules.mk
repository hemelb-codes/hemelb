include $(MK)/header.mk

SRCS := HttpPost.cc \
	Network.cc \
	SimulationParameters.cc \
	steering.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
