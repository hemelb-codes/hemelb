include $(MK)/header.mk

SRCS := BaseTopologyManager.cc

ifdef HEMELB_CFG_MULTIMACHINE
  SRCS += MultiMachineTopologyManager.cc
else
  SRCS += SingleMachineTopologyManager.cc
endif

TARGETS := libHemeLbTopology.$(LIBEXT)

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
