#CC=mpicxx
#CC=g++
#CC=mpicxx

CC=lamc++-64

INCLUDES=-I/usr/local/include
#LIBS=-L/usr/local/lib -llammpio -llammpi++ -lmpi -llam -ldl
EXE=hemelb

MPIRUN=/usr/local/bin/mpiexec
MPIRUN=/usr/local/bin/mpirun
MPIRUN=/usr/bin/mpirun
MPIRUN=/usr/bin/mpirun

MPICPUS=2

TEST1INPUT=/Users/stevenmanos/hemelb/Input/square_duct_32x16x16_input.asc
TEST2INPUT=/Users/stevenmanos/hemelb/Input/angio1_input.asc

TEST1INPUT=/home/steven/hemelb/Input/square_duct_32x16x16_input.asc
TEST2INPUT=/home/marco/hemelb/Input/angio1_input.asc #/home/steven/hemelb/Input/angio1_input.asc

FLAGS=-O3 ${REG_CPPFLAGS}
DEFS=-DRG #-DSTEER

OBJ=config.o \
	io.o \
	lb.o \
	main.o \
	rt.o \
	topology.o \
	eVizRLEUtil.o
	
$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(LIBS) ${REG_LDFLAGS} $(OBJ) -o $(EXE)

test1:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST1INPUT)

test2:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST2INPUT)

display:
#	g++ -O2 -framework OpenGL -framework GLUT -framework Foundation -o displayClient displayClient.cc eVizRLEUtil.c 
	g++ displayClient.cc -O3 -lGLU -lglut eVizRLEUtil.c -o displayClient

.cc.o:
	$(CC) $(FLAGS) $(DEFS) $(INCLUDES) -c $<

.c.o:
	$(CC) $(FLAGS) $(DEFS) -c $<

clean:
	rm -f $(OBJ) $(EXE) displayClient

distclean: clean
	find ./ -name \*~ -exec rm -f {} \;

