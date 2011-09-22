include $(MK)/header.mk

SRCS := BlockTraverser.cc \
	BlockTraverserWithVisitedBlockTracker.cc \
	GeometryReader.cc \
        GlobalLatticeData.cc \
        LatticeData.cc \
	LocalLatticeData.cc \
	SiteTraverser.cc \
	VolumeTraverser.cc
        

TARGETS := libHemeLbGeometry.$(LIBEXT)

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
