CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		#-DNDEBUG
OBJS=		utils.o seq.o sais.o saux.o rld.o exact.o cmd.o
PROG=		fermi
INCLUDES=	
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

fermi:$(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

rld.o:rld.h
cmd.o:fermi.h rld.h
exact.o:fermi.h rld.h kstring.h

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM
