include $(MK)/header.mk

SUBDIRS := xdr ascii

TARGETS = Writer.o 
	  
INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk