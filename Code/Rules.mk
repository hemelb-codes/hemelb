# -*- mode: makefile;-*-
include $(MK)/header.mk

TARGETS = $(EXE) hlbTest
SUBDIRS = vis steering lb io debug topology xml parmetis

$(EXE)_DEPS = D3Q15.o \
	fileutils.o \
	net.o \
        SimulationMaster.o \
        SimConfig.o \
	usage.o \
        utilityFunctions.o \
	main.o \
	$(SUBDIRS_TGTS)

hlbTest_DEPS = D3Q15.o \
        fileutils.o \
        net.o \
        SimulationMaster.o \
        SimConfig.o \
        usage.o \
        utilityFunctions.o \
        tests.o \
        $(SUBDIRS_TGTS)

HEMELB_INCLUDEPATHS += $(TOP)

include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
