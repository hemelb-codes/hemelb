include $(MK)/header.mk

TARGETS := libHemeLbSteering.$(LIBEXT)

SRCS :=

ifeq ($(HEMELB_STEERING_LIB), none)
# no steering - use "off"
SUBDIRS := none

else ifeq ($(HEMELB_STEERING_LIB), basic)
# Steering enabled, use "on"
SUBDIRS := basic

endif

SUBDIRS += common

HEMELB_DEFS += HEMELB_STEERING_LIB=$(HEMELB_STEERING_LIB) HEMELB_STEERING_LIB_$(HEMELB_STEERING_LIB)

$(TARGETS)_DEPS = $(SUBDIRS_TGTS)

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
