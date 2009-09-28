include make/Makefile.${HEMELB_MACHINE}

OBJ=config.o \
	io.o \
	lb.o \
	main.o \
	rt.o \
	topology.o \
	network.o \
	benchmark.o \
	colourpalette.o \
	visthread.o \
	fileutils.o \
	usage.o \
	steering-common.o \
	steering.o \
	steering-sim-params.o \
	http_post.o 
#	realitygrid.o

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) $(LIBS) -o $(EXE)

doc:
	doxygen 

display:
	g++ displayClient.cc -O3 -lGLU -lglut eVizRLEUtil.c -o displayClient

mac-display:
	g++ -O3 -framework OpenGL -framework GLUT -framework Foundation -o displayClient displayClient.cc eVizRLEUtil.c -DMAC_GLUT

display-debug:
	g++ displayClient.cc -g -lGLU -lglut eVizRLEUtil.c -o displayClient

.cc.o:
	$(CC) $(FLAGS) $(DEFS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJ) $(EXE) displayClient env_details.asc

distclean: clean
	find ./ -name \*~ -exec rm -f {} \;

include make/Makefile.Tests

