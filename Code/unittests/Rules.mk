include $(MK)/header.mk
ifndef NOTESTS
TARGETS := libHemeLbUnitTests.$(LIBEXT)

$(TARGETS)_DEPS = main.o

INCLUDES_$(d) := $(INCLUDES_$(parent))
endif
include $(MK)/footer.mk
