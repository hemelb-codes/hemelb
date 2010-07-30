include make/Makefile.${HEMELB_MACHINE}

OBJ=config.o \
        xdrReader.o \
        xdrWriter.o \
        xdrFileWriter.o \
        xdrMemWriter.o \
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
#	realitygrid.o

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) $(LIBS) -o $(EXE)

doc:
	doxygen 

.cc.o:
	$(CC) $(FLAGS) $(DEFS) $(INCLUDES) -c $<

clean:
	-rm -f $(OBJ) $(EXE) env_details.asc
	cd client && make clean

distclean: clean
	find ./ -name \*~ -exec rm -f {} \;

include make/Makefile.Tests

