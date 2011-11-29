include $(MK)/header.mk

TARGETS = AsciiFileWriter.o AsciiStreamWriter.o
	  
INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
