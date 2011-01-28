include $(MK)/header.mk

TARGETS := OsxDebugger.o

INCLUDES_$(d) := $(INCLUDES_$(parent))
include $(MK)/footer.mk
