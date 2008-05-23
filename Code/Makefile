include make/Makefile.${MACHINE}

OBJ=config.o \
	io.o \
	lb.o \
	main.o \
	rt.o \
	topology.o 

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(LIBS) ${REG_LDFLAGS} $(OBJ) -o $(EXE)

.cc.o:
	$(CC) $(FLAGS) $(DEFS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJ) $(EXE) displayClient

distclean: clean
	find ./ -name \*~ -exec rm -f {} \;

