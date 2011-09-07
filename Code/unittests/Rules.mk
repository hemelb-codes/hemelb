include $(MK)/header.mk

TARGETS := libHemeLbUnitTests.$(LIBEXT)

SUBDIRS := lbtests

$(TARGETS)_DEPS = $(SUBDIRS_TGTS) \
                  main.o \
                  $(foreach sd,$(SUBDIRS_$(d)),$(call subtree_tgts,$(sd)))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
