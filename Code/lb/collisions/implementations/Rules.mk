include $(MK)/header.mk

TARGETS := Implementation.o \

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
