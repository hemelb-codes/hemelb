CC=lamc++-64
EXE=hemelb
MPIRUN=/usr/bin/mpirun
MPICPUS=2
TEST1INPUT=/home/steven/hemelb/Input/square_duct_32x16x16_input.asc
TEST2INPUT=/home/steven/hemelb/Input/angio1_input.asc

FLAGS=-O3
DEFS=-DRG

OBJ=config.o \
	io.o \
	lb.o \
	main.o \
	rt.o \
	topology.o

$(EXE):$(OBJ)
	$(CC) $(FLAGS) $(OBJ) -o $(EXE)

test1:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST1INPUT)

test2:
	$(MPIRUN) -np $(MPICPUS) ./$(EXE) $(TEST2INPUT)

.cc.o:
	$(CC) $(FLAGS) $(DEFS) -c $<

clean:
	rm -f $(OBJ) $(EXE)

