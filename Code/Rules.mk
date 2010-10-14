# -*- mode: makefile;-*-
include $(MK)/header.mk

TARGETS = $(EXE)
SUBDIRS = vis steering io dbg lb

$(EXE)_DEPS = config.o \
        D3Q15.o \
	io.o \
	lb.o \
	fileutils.o \
	net.o \
	usage.o \
        utilityFunctions.o \
	main.o \
	$(SUBDIRS_TGTS)


HEMELB_INCLUDEPATHS += $(TOP)

include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
