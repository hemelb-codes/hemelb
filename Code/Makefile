CC=mpicxx
EXE=hemelb
MPIRUN=/usr/local/bin/mpiexec
MPICPUS=2
TEST1INPUT=/home/steven/hemelb/Input/square_duct_32x16x16_input.asc
TEST2INPUT=/home/steven/hemelb/Input/angio1_input.asc

REG_HOME=/home/ngs0011/work/RealityGrid/file/reg_steer_lib
REG_CPPFLAGS=-I${REG_HOME}/include
REG_LDFLAGS=-L${REG_HOME}/lib64 -Wl,-Bstatic -lReG_Steer -lReG_Steer_Utils \
		-lReG_Steer_SOAP -Wl,-Bdynamic -lxml2

FLAGS=-O3 ${REG_CPPFLAGS}
DEFS=-DRG -DSTEER

OBJ=config.o \
	io.o \
	lb.o \
	main.o \
	rt.o \
	topology.o \
	eVizRLEUtil.o

$(EXE): $(OBJ)
	$(CC) $(FLAGS) ${REG_LDFLAGS} $(OBJ) -o $(EXE)

test1:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST1INPUT)

test2:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST2INPUT)

display:
	g++ displayClient.cc -O2 -lGLU -lglut eVizRLEUtil.c -o displayClient

.cc.o:
	$(CC) $(FLAGS) $(DEFS) -c $<

.c.o:
	$(CC) $(FLAGS) $(DEFS) -c $<

clean:
	rm -f $(OBJ) $(EXE) displayClient

distclean: clean
	find ./ -name \*~ -exec rm -f {} \;
