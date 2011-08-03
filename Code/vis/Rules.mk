include $(MK)/header.mk

TARGETS := libHemeLbVis.$(LIBEXT)

SUBDIRS := rayTracer

$(TARGETS)_DEPS = $(SUBDIRS_TGTS) \
	ColPixel.o \
        GlyphDrawer.o \
	StreaklineDrawer.o \
        Control.o \
	Screen.o \
        ScreenPixels.o \
        Viewpoint.o 

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
