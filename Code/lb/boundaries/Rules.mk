include $(MK)/header.mk

TARGETS = BoundaryValues.o \
	  BoundaryComms.o \
	  InOutLet.o \
	  InOutLetCosine.o \
	  InOutLetFile.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
