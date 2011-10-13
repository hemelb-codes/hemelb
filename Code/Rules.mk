# -*- mode: makefile;-*-
include $(MK)/header.mk

TARGETS = $(EXE) $(UNITTESTS)

# Note that ParMetis doesn't count as a subdirectory because it uses its own
# build system.
# To rebuild that, go to the ParMetis directory and run make there.

SUBDIRS = steering vis lb net debug topology xml util geometry io log

$(EXE)_DEPS = D3Q15.o \
        SimulationMaster.o \
        SimConfig.o \
        main.o \
        mpiInclude.o \
        $(SUBDIRS_TGTS)

SUBDIRS = unittests steering vis lb net debug topology xml util geometry io log

$(UNITTESTS)_DEPS = D3Q15.o \
                    SimulationMaster.o \
                    SimConfig.o \
                    mpiInclude.o \
                    $(SUBDIRS_TGTS)

include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
