include $(MK)/header.mk

SRCS := LoggerCommon.cc


ifeq ($(HEMELB_LOG_LEVEL), debug)
    SRCS += LoggerDebug.cc
else ifeq ($(HEMELB_LOG_LEVEL), warning)
    SRCS += LoggerWarning.cc
else
    SRCS += LoggerInfo.cc
endif

TARGETS := libHemeLbLogger.$(LIBEXT)

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
