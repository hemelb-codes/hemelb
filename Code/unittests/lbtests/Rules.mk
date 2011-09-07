include $(MK)/header.mk

TARGETS := libHemeLbUnitTestsLb.$(LIBEXT)

$(TARGETS)_DEPS = KernelTests.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
