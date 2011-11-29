include $(MK)/header.mk

TARGETS := libHemeLbIo.$(LIBEXT)

SUBDIRS := writers

$(TARGETS)_DEPS = $(foreach sd,$(SUBDIRS_$(d)),$(call subtree_tgts,$(sd)))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
