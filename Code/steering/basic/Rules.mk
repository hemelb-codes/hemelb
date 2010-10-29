include $(MK)/header.mk

SRCS := HttpPost.cc \
	Network.cc \
	SimulationParameters.cc \
	steering.cc \
	Control.cc \
	Threadable.cc \
	NetworkThread.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
