include $(MK)/header.mk

SRCS := writer.cc \
	xdrReader.cc \
	xdrWriter.cc \
	xdrFileWriter.cc \
	xdrMemWriter.cc \
	asciiStreamWriter.cc \
	asciiFileWriter.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
