include $(MK)/header.mk

#SRCS := NullDebugger.cc
TARGETS := NullDebugger.o

INCLUDES_$(d) := $(INCLUDES_$(parent))
include $(MK)/footer.mk
