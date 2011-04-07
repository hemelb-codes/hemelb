include $(MK)/header.mk

SRCS := LocalLatticeData.cc \
        GlobalLatticeData.cc \
        LatticeData.cc \
        GeometryReader.cc

TARGETS := libHemeLbGeometry.$(LIBEXT)

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
