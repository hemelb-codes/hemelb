include $(MK)/header.mk

include $(MK)/Makefile.${HEMELB_MACHINE}

TARGETS = $(EXE)
SUBDIRS = io

hemelb_DEPS = config.o \
	io.o \
	lb.o \
        visualisationControl.o \
        glyphDrawer.o \
        streaklineDrawer.o \
        rayTracer.o \
	benchmark.o \
	http_post.o \
	fileutils.o \
        colpixel.o \
        rt.o \
	net.o \
	network.o \
	colourpalette.o \
	visthread.o \
	usage.o \
	steering-common.o \
	steering.o \
	steering-sim-params.o \
        utilityFunctions.o \
	main.o \
	$(SUBDIRS_TGTS)
app.exe_LIBS = $(LIBS)


include $(MK)/footer.mk
# This is just a convenience - to let you know when make has stopped
# interpreting make files and started their execution.
$(info Rules generated...)

# vim: set ft=make :
