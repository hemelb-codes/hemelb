# -*- mode: makefile;-*-
include $(MK)/header.mk

TARGETS = $(EXE)

# Note that ParMetis doesn't count as a subdirectory because it uses its own
# build system.
# To rebuild that, go to the ParMetis directory and run make there.

SUBDIRS = geometry steering vis lb net debug topology io xml util

$(EXE)_DEPS = D3Q15.o \
        SimulationMaster.o \
        SimConfig.o \
        main.o \
        $(SUBDIRS_TGTS)

include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
