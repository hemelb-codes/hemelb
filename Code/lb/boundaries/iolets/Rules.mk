include $(MK)/header.mk

TARGETS = InOutLet.o \
	  InOutLetCosine.o \
	  InOutLetFile.o
	  
INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
