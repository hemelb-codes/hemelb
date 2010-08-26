include $(MK)/header.mk

TARGETS := libHemeLbIo.$(LIBEXT)
SRCS := Writer.cc \
	XdrReader.cc \
	XdrWriter.cc \
	XdrFileWriter.cc \
	XdrMemWriter.cc \
	AsciiStreamWriter.cc \
	AsciiFileWriter.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
