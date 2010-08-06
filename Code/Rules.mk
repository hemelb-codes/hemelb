# -*- mode: makefile;-*-
include $(MK)/header.mk

TARGETS = $(EXE)
SUBDIRS = io vis

$(EXE)_DEPS = config.o \
	io.o \
	lb.o \
	benchmark.o \
	http_post.o \
	fileutils.o \
	net.o \
	network.o \
	usage.o \
	steering-common.o \
	steering.o \
	steering-sim-params.o \
        utilityFunctions.o \
	main.o \
	$(SUBDIRS_TGTS)

HEMELB_INCLUDEPATHS += $(TOP)

include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
