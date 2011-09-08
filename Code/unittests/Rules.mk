include $(MK)/header.mk

TARGETS := libHemeLbUnitTests.$(LIBEXT)

$(TARGETS)_DEPS = main.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
