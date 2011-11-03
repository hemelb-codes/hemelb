include $(MK)/header.mk

TARGETS := libConfiguration.$(LIBEXT)
SRCS := SimConfig.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
