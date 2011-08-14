CC=			gcc
CFLAGS=		-g -Wall #-O2
DFLAGS=		#-D_USE_RLE6 #-DNDEBUG
OBJS=		utils.o seq.o sais.o sais64.o saux.o rld.o bit3.o exact.o merge.o \
			append.o build.o cmd.o
PROG=		fermi
INCLUDES=	
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

fermi:$(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

sais64.o:sais.c
		$(CC) -c $(CFLAGS) -D_SAIS64 $(DFLAGS) sais.c -o $@

bit3.o:bit3.h
rld.o:rld.h bit3.h
build.o:fermi.h rld.h
append.o:fermi.h rld.h
exact.o:fermi.h rld.h kstring.h kvec.h
merge.o:fermi.h rld.h khash.h ksort.h
cmd.o:fermi.h rld.h bit3.h kseq.h
main.o:fermi.h

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session*
