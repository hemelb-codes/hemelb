include $(MK)/header.mk

TARGETS := libHemeLbVis.$(LIBEXT)

SUBDIRS := rayTracer streaklineDrawer

$(TARGETS)_DEPS = $(SUBDIRS_TGTS) \
                  GlyphDrawer.o \
                  Control.o \
                  Screen.o \
                  Viewpoint.o \
                  BasicPixel.o \
                  Rendering.o \
                  ResultPixel.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
