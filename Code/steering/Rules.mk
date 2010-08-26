include $(MK)/header.mk

TARGETS := libHemeLbSteering.$(LIBEXT)

SRCS := common.cc

ifeq ($(HEMELB_STEERING_LIB), none)
# no steering - use "off"
SUBDIRS := none

else ifeq ($(HEMELB_STEERING_LIB), basic)
# Steering enabled, use "on"
SUBDIRS := basic

endif

SUBDIRS += common

$(TARGETS)_DEPS = $(SUBDIRS_TGTS)

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
