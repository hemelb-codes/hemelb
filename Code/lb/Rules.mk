include $(MK)/header.mk

TARGETS := libHemeLbLb.$(LIBEXT)

$(TARGETS)_DEPS = lb.o \
                  io.o \
                  HFunction.o \
                  StabilityTester.o \
                  EntropyTester.o \
                  SimulationState.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
