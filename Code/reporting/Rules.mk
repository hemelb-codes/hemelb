include $(MK)/header.mk

TARGETS := libHemeLbReporting.$(LIBEXT)
SRCS := FileManager.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
