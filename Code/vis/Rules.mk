include $(MK)/header.mk

TARGETS := libHemeLbVis.$(LIBEXT)

SUBDIRS := rayTracer streaklineDrawer

$(TARGETS)_DEPS = $(SUBDIRS_TGTS) \
	BlockTraverser.o \
	ColPixelMPIDataTypes.o \
	GlyphDrawer.o \
	Control.o \
	Screen.o \
	SiteTraverser.o \
	Viewpoint.o \
	VolumeTraverser.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
