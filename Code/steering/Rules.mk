include $(MK)/header.mk

SRCS := SimulationParameters.cc \
	common.cc \
	HttpPost.cc \
	Network.cc \
	steering.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
