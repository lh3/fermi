CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		-D_DNA_ONLY #-D_USE_RLE6 #-DNDEBUG
OBJS=		utils.o seq.o ksa.o ksa64.o rld.o exact.o merge.o correct.o \
			build.o smem.o unitig.o msg.o bog.o cmd.o
PROG=		fermi
INCLUDES=	
LIBS=		-lpthread -lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

fermi:$(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

ksa64.o:ksa.c
		$(CC) -c $(CFLAGS) -D_KSA64 $(DFLAGS) ksa.c -o $@

rld.o:rld.h
build.o:fermi.h rld.h
exact.o:fermi.h rld.h kstring.h kvec.h
unitig.o:fermi.h rld.h kstring.h kvec.h
correct.o:fermi.h rld.h kvec.h
smem.o:fermi.h rld.h kvec.h
merge.o:fermi.h rld.h ksort.h
bog.o:fermi.h kstring.h
cmd.o:fermi.h rld.h kseq.h
main.o:fermi.h

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*
