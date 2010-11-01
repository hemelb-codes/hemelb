include $(MK)/header.mk

TARGETS := libXml.$(LIBEXT)
SRCS := action_store.cc \
        htmlutil.cc \
        lex_util.cc \
        node_set.cc \
        ticpp.cc \
        tinystr.cc \
        tinyxml.cc \
        tinyxmlerror.cc \
        tinyxmlparser.cc \
        tokenlist.cc \
        xml_util.cc \
        xpath_expression.cc \
        xpath_processor.cc \
        xpath_stack.cc \
        xpath_static.cc \
        xpath_stream.cc \
        xpath_syntax.cc

$(TARGETS)_DEPS := $(subst .cc,.$(OBJEXT), $(SRCS))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
