include $(MK)/header.mk

TARGETS := libHemeLbDbg.$(LIBEXT)
#SRCS := Debugger.cc

ifeq ($(HEMELB_DBG_LEVEL), 0)
# I.e. no debugging
SUBDIRS := none
HEMELB_DBG_LIBRARY := none

else # HEMELB_DBG_LEVEL nonzero
$(info Building in debug mode.)

ifdef HEMELB_CFG_ON_OSX
SUBDIRS := common OSX
HEMELB_DBG_LIBRARY := OSX

else ifdef HEMELB_CFG_ON_LINUX
SUBDIRS := common linux
HEMELB_DBG_LIBRARY := linux

else # non OSX/Linux debug interfaces not implemented
SUBDIRS := none
HEMELB_DBG_LIBRARY := none

endif

endif # HEMELB_DBG_LEVEL==0

HEMELB_DEFS += HEMELB_DBG_LIBRARY=$(HEMELB_DBG_LIBRARY) HEMELB_DBG_LIBRARY_$(HEMELB_DBG_LIBRARY)

$(TARGETS)_DEPS = Debugger.o $(SUBDIRS_TGTS)

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
