include $(MK)/header.mk

SRCS := Writer.cc \
	XdrReader.cc \
	XdrWriter.cc \
	XdrFileWriter.cc \
	XdrMemWriter.cc \
	AsciiStreamWriter.cc \
	AsciiFileWriter.cc

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
