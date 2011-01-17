include $(MK)/header.mk

TARGETS := libHemeLbUtil.$(LIBEXT)
SRCS := fileutils.cc \
        utilityFunctions.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
