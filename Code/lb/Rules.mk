include $(MK)/header.mk

TARGETS := libHemeLbLb.$(LIBEXT)

SUBDIRS := collisions streamers rheology_models

$(TARGETS)_DEPS = lb.o \
                  io.o \
                  HFunction.o \
                  StabilityTester.o \
                  EntropyTester.o \
                  SimulationState.o \
                  $(foreach sd,$(SUBDIRS_$(d)),$(call subtree_tgts,$(sd)))

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
