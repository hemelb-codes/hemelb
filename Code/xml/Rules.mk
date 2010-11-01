include $(MK)/header.mk

TARGETS := libXml.$(LIBEXT)
SRCS := ticpp.cc \
        tinystr.cc \
        tinyxml.cc \
        tinyxmlerror.cc \
        tinyxmlparser.cc \

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
