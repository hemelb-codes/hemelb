include $(MK)/header.mk

SRCS := colPixel.cc \
	colourPalette.cc \
	glyphDrawer.cc \
	rayTracer.cc \
	rt.cc \
	streaklineDrawer.cc \
	visthread.cc \
	visualisationControl.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
