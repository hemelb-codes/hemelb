include ../Makefile.in

BINDIR = ../Graphs

INCLUDES = -I./ -I../ParMETISLib $(INCDIR) 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)


LIBSDIR = -L.. $(LIBDIR) 
LIBS = -lparmetis -lmetis $(XTRALIBS) -lm
PARMETISLIB = ../libparmetis.a ../libmetis.a


PTESTOBJS = ptest.o io.o adaptgraph.o 
MESHTESTOBJS = mtest.o io.o
PARMETISOBJS = parmetis.o io.o adaptgraph.o
POMETISOBJS = pometis.o io.o 

          
.c.o:
	$(CC) $(CFLAGS) -c $*.c


default: $(BINDIR)/ptest$(VERNUM) $(BINDIR)/mtest$(VERNUM) $(BINDIR)/parmetis$(VERNUM) \
         $(BINDIR)/pometis$(VERNUM)

$(BINDIR)/mtest$(VERNUM): $(MESHTESTOBJS) $(PARMETISLIB)
	$(LD)  -o $@ $(MESHTESTOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@

$(BINDIR)/ptest$(VERNUM): $(PTESTOBJS) $(PARMETISLIB)
	$(LD)  -o $@ $(PTESTOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@

$(BINDIR)/parmetis$(VERNUM): $(PARMETISOBJS) $(PARMETISLIB)
	$(LD)  -o $@ $(PARMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@

$(BINDIR)/pometis$(VERNUM): $(POMETISOBJS) $(PARMETISLIB)
	$(LD)  -o $@ $(POMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@


clean:
	rm -f *.o ;\
	rm -f $(BINDIR)/ptest$(VERNUM) 
	rm -f $(BINDIR)/mtest$(VERNUM) 
	rm -f $(BINDIR)/parmetis$(VERNUM) 
	rm -f $(BINDIR)/pometis$(VERNUM) 

realclean:
	rm -f *.o ;\
	rm -f $(BINDIR)/ptest$(VERNUM) 
	rm -f $(BINDIR)/mtest$(VERNUM) 
	rm -f $(BINDIR)/parmetis$(VERNUM) 
	rm -f $(BINDIR)/pometis$(VERNUM) 

