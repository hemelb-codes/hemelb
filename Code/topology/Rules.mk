include $(MK)/header.mk

SRCS := NetworkTopology.cc

ifdef HEMELB_CFG_MULTIMACHINE
  SRCS += MultiMachineNetworkTopology.cc
else
  SRCS += SingleMachineNetworkTopology.cc
endif

TARGETS := libHemeLbTopology.$(LIBEXT)

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
