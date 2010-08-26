include $(MK)/header.mk

SRCS := debug.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))
include $(MK)/footer.mk
