CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=
OBJS=		sais.o rld.o index.o
PROG=		bwa2
INCLUDES=	
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwa2:$(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

rld.o:rld.h

cleanlocal:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

clean:cleanlocal-recur
