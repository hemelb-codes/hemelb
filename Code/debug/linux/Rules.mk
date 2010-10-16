include $(MK)/header.mk

TARGETS := LinuxDebugger.o

INCLUDES_$(d) := $(INCLUDES_$(parent))
include $(MK)/footer.mk
