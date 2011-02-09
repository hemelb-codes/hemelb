include $(MK)/header.mk

TARGETS := libHemeLbMethod.$(LIBEXT)

#SRCS := lb.cc \
         io.cc \
         StabilityTester.cc         

SUBDIRS := collisions

$(TARGETS)_DEPS = $(SUBDIRS_TGTS) \
                   lb.o \
                   io.o \
                   LocalLatticeData.o

INCLUDES_$(d) := $(INCLUDES_$(parent))

include $(MK)/footer.mk
