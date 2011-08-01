CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		#-DNDEBUG
OBJS=		sais.o saux.o rle6.o rld.o index.o
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
rle6.o:rle6.h
index.o:rle6.h rld.h

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM
