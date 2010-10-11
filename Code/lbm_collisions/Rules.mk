include $(MK)/header.mk

TARGETS := libHemeLbCollisions.$(LIBEXT)
SRCS := Collisions.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
