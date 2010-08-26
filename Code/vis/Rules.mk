include $(MK)/header.mk

TARGETS := libHemeLbVis.$(LIBEXT)
SRCS := ColPixel.cc \
	ColourPalette.cc \
	Layer.cc \
	GlyphDrawer.cc \
	RayTracer.cc \
	StreaklineDrawer.cc \
	visthread.cc \
	Control.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
