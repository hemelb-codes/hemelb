# -*- mode: makefile;-*-
include $(MK)/header.mk

TARGETS = $(EXE) hlbTest
SUBDIRS = vis steering lb io debug topology xml parmetis util

$(EXE)_DEPS = D3Q15.o \
	net.o \
        SimulationMaster.o \
        SimConfig.o \
	main.o \
	$(SUBDIRS_TGTS)

hlbTest_DEPS = D3Q15.o \
        net.o \
        SimulationMaster.o \
        SimConfig.o \
        tests.o \
        $(SUBDIRS_TGTS)

HEMELB_INCLUDEPATHS += $(TOP)

include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
