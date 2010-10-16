include $(MK)/header.mk

#SRCS := NullDebugger.cc
TARGETS := ActiveDebugger.o

INCLUDES_$(d) := $(INCLUDES_$(parent))
include $(MK)/footer.mk
