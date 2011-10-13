include $(MK)/header.mk

TARGETS := libHemeLbVis.$(LIBEXT)

SUBDIRS := rayTracer

$(TARGETS)_DEPS = $(SUBDIRS_TGTS) \
                  GlyphDrawer.o \
                  StreaklineDrawer.o \
                  Control.o \
                  Screen.o \
                  Viewpoint.o \
                  BasicPixel.o \
                  StreakPixel.o \
                  Rendering.o \
                  ResultPixel.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
