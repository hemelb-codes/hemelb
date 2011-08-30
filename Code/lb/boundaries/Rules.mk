include $(MK)/header.mk

TARGETS = BoundaryValues.o \
	  BoundaryComms.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
