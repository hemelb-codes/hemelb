include $(MK)/header.mk

SRCS := TopologyManager.cc \

TARGETS := libHemeLbTopology.$(LIBEXT)

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
