CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		#-D_USE_RLE6 #-DNDEBUG
OBJS=		utils.o seq.o ksa.o ksa64.o rld.o exact.o merge.o sub.o correct.o \
			build.o smem.o unitig.o seqsort.o cmp.o cmd.o example.o \
			ksw.o mag.o bubble.o scaf.o bcr.o bprope6.o ropebwt.o
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

rld.o:rld.c rld.h
		$(CC) -c $(CFLAGS) -D_DNA_ONLY $(DFLAGS) $(INCLUDES) $< -o $@

build.o:build.c fermi.h rld.h
exact.o:exact.c fermi.h rld.h kstring.h kvec.h
unitig.o:unitig.c fermi.h rld.h kstring.h kvec.h
correct.o:correct.c fermi.h rld.h kvec.h kseq.h
smem.o:smem.c fermi.h rld.h kvec.h
merge.o:merge.c fermi.h rld.h ksort.h
sub.o:sub.c fermi.h rld.h
cmd.o:cmd.c fermi.h rld.h kseq.h
mag.o:mag.c mag.h kseq.h
bubble.o:bubble.c mag.h ksw.h
scaf.o:scaf.c mag.h rld.h fermi.h kvec.h khash.h ksw.h
cmp.o:cmp.c rld.h fermi.h kvec.h
main.o:main.c fermi.h

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*
