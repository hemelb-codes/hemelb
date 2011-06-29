include $(MK)/header.mk

TARGETS := libXml.$(LIBEXT)
SRCS := tinystr.cpp \
        tinyxml.cpp \
        tinyxmlerror.cpp \
        tinyxmlparser.cpp

$(TARGETS)_DEPS := $(subst .cpp,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
