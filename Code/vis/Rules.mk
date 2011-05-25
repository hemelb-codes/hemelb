include $(MK)/header.mk

TARGETS := libHemeLbVis.$(LIBEXT)
SRCS := ColPixel.cc \
        GlyphDrawer.cc \
        RayTracer.cc \
        StreaklineDrawer.cc \
        Control.cc \
        Screen.cc \
        ScreenPixels.cc \
        Viewpoint.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
