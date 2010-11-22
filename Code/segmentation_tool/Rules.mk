# -*- mode: makefile;-*-
include $(MK)/header.mk

TARGETS = $(EXE)
SUBDIRS =

$(EXE)_DEPS = config.o \
              editing.o \
              io.o \
              main.o \
	      io_dicom.o \
	      math.o \
              menu.o \
              rt.o \
              segmentation.o \
              timing.o \
   	      usage.o \
              vis.o

$(EXE)_INCLUDEPATHS += $(TOP)

include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
