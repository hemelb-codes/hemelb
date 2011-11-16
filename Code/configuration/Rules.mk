include $(MK)/header.mk

TARGETS := libHemeLbConfiguration.$(LIBEXT)
SRCS := SimConfig.cc \
		CommandLine.cc 

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
