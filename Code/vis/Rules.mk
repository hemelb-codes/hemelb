include $(MK)/header.mk

SRCS := ColPixel.cc \
	ColourPalette.cc \
	GlyphDrawer.cc \
	RayTracer.cc \
	StreaklineDrawer.cc \
	visthread.cc \
	Control.cc


INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
