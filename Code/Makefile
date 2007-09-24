#CC=mpicxx
#CC=g++

CC=lamc++-64

#INCLUDES=-I/usr/local/include
#LIBS=-L/usr/local/lib -llammpio -llammpi++ -lmpi -llam -ldl
EXE=hemelb

MPIRUN=/usr/local/bin/mpiexec
MPIRUN=/usr/local/bin/mpirun
MPIRUN=/usr/bin/mpirun
MPIRUN=/usr/bin/mpirun

#REG_HOME=/home/steven/reg_steer_lib

REG_CPPFLAGS=-I${REG_HOME}/include

#REG_LDFLAGS=-L${REG_HOME}/lib64 -Wl,-Bstatic -lReG_Steer -lReG_Steer_Utils \
#                -lReG_Steer_SOAP -Wl,-Bdynamic -lxml2

MPICPUS=2

TEST1INPUT=$(HOME)/hemelb/Input/square_duct_32x16x16
TEST2INPUT=$(HOME)/hemelb/Input/angio1
TEST3INPUT=$(HOME)/hemelb/Input/bifurcation_512x512x512

FLAGS=-O3 ${REG_CPPFLAGS} 
#DEFS=-DRG -DSTEER
DEFS=-DRG 

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

test3:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST3INPUT)

display:
	g++ displayClient.cc -O3 -lGLU -lglut eVizRLEUtil.c -o displayClient

mac-display:
	g++ -O3 -framework OpenGL -framework GLUT -framework Foundation -o displayClient displayClient.cc eVizRLEUtil.c 

display-debug:
	g++ displayClient.cc -g -lGLU -lglut eVizRLEUtil.c -o displayClient

.cc.o:
	$(CC) $(FLAGS) $(DEFS) $(INCLUDES) -c $<

.c.o:
	$(CC) $(FLAGS) $(DEFS) -c $<

clean:
	rm -f $(OBJ) $(EXE) displayClient

distclean: clean
	find ./ -name \*~ -exec rm -f {} \;

