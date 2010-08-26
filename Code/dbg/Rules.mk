include $(MK)/header.mk

TARGETS := libHemeLbDbg.$(LIBEXT)

ifeq ($(HEMELB_DBG_LEVEL), 0)
# I.e. no debugging
SUBDIRS := none
# SRCS := none/debug.cc

else # HEMELB_DBG_LEVEL nonzero

ifdef $(HEMELB_CFG_ON_OSX)
SUBDIRS := OSX
# SRCS := OSX/debug.cc

else # non OSX debug interfaces not implemented
SUBDIRS := none
endif

endif # HEMELB_DBG_LEVEL==0


$(TARGETS)_DEPS = $(SUBDIRS_TGTS)

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
