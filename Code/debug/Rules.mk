include $(MK)/header.mk

TARGETS := libHemeLbDebug.$(LIBEXT)
#SRCS := Debugger.cc

ifeq ($(HEMELB_DEBUG_LEVEL), 0)
# I.e. no debugging
SUBDIRS := none
HEMELB_DEBUG_LIBRARY := none

else # HEMELB_DEBUG_LEVEL nonzero
$(info Building in debug mode.)

ifdef HEMELB_CFG_ON_OSX
SUBDIRS := common OSX
HEMELB_DEBUG_LIBRARY := OSX

else ifdef HEMELB_CFG_ON_LINUX
SUBDIRS := common linux
HEMELB_DEBUG_LIBRARY := linux

else # non OSX/Linux debug interfaces not implemented
SUBDIRS := none
HEMELB_DEBUG_LIBRARY := none

endif

endif # HEMELB_DEBUG_LEVEL==0

HEMELB_DEFS += HEMELB_DEBUG_LIBRARY=$(HEMELB_DEBUG_LIBRARY) HEMELB_DEBUG_LIBRARY_$(HEMELB_DEBUG_LIBRARY)

$(TARGETS)_DEPS = Debugger.o $(SUBDIRS_TGTS)

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
