include $(MK)/header.mk

TARGETS = XdrFileReader.o XdrFileWriter.o XdrMemReader.o XdrMemWriter.o XdrWriter.o XdrReader.o
	  
INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
