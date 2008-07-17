#/home/steven/lammpi/bin/lamboot -v ~/lammpinodes
#/home/steven/lammpi/bin/lamhalt

CC=icpc

EXE=hemelb
#LIBS=-L/usr/lib64/lam -llammpio -llammpi++ -llamf77mpi -lmpi -llam -laio -laio -lutil -ldl 

MPIRUN=/home/steven/lammpi/bin/mpirun
INCLUDES=-I/home/steven/lammpi/include

#REG_HOME=/home/steven/reg_steer_lib

REG_CPPFLAGS=-I${REG_HOME}/include

#REG_LDFLAGS=-L${REG_HOME}/lib64 -Wl,-Bstatic -lReG_Steer -lReG_Steer_Utils \
#                -lReG_Steer_SOAP -Wl,-Bdynamic -lxml2  1000  1000  0.001

MPICPUS=4

TEST1INPUT=$(HOME)/hemelb/Input/square_duct_32x16x16     2 10000000 0.01
TEST2INPUT=$(HOME)/hemelb/Input/angio1                   10 200000 0.0001171875
TEST3INPUT=$(HOME)/hemelb/Input/bifurcation_512x512x512  10 200000 0.0001
TEST4INPUT=$(HOME)/hemelb/Input/angio4                   1 2000 0.0009375
TEST5INPUT=$(HOME)/hemelb/Input/angio2                   10 200000 0.00234375
TEST6INPUT=$(HOME)/hemelb/Input/angio8                   10 200000 0.0002245

FLAGS=-O3 -xT -ipo -pthread -Wunused
FLAGS=-pthread
DEFS=

LIBS=-L/home/steven/lammpi/lib -llammpio -llammpi++ -lmpi -llam -ldl

OBJ=config.o \
	io.o \
	lb.o \
	main.o \
	rt.o \
	topology.o \
	steering.o \
	network.o \
	benchmark.o \
	colourpalette.o \
	visthread.o \
	usage.o

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) $(LIBS) -o $(EXE)

test1:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST1INPUT)

test2:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST2INPUT)

test3:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST3INPUT)

test4:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST4INPUT)

test5:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST5INPUT)

test6:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST6INPUT)

display:
	g++ displayClient.cc -O3 -lGLU -lglut eVizRLEUtil.c -o displayClient

mac-display:
	g++ -O3 -framework OpenGL -framework GLUT -framework Foundation -o displayClient displayClient.cc eVizRLEUtil.c -DMAC_GLUT

display-debug:
	g++ displayClient.cc -g -lGLU -lglut eVizRLEUtil.c -o displayClient

.cc.o:
	$(CC) $(FLAGS) $(DEFS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJ) $(EXE) displayClient

distclean: clean
	find ./ -name \*~ -exec rm -f {} \;
