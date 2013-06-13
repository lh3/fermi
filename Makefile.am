bin_PROGRAMS = fermi
fermi_SOURCES = utils.c seq.c ksa.c exact.c merge.c sub.c correct.c \
                build.c smem.c unitig.c seqsort.c cmp.c cmd.c example.c \
                ksw.c mag.c bubble.c scaf.c bcr.c bprope6.c ropebwt.c main.c
fermi_LDADD = libksa64.a librld.a
noinst_LIBRARIES = libksa64.a librld.a
libksa64_a_SOURCES = ksa.c
libksa64_a_CFLAGS ="-D_KSA64 $(CFLAGS)"
librld_a_SOURCES = rld.c rld.h
librld_a_CFLAGS = "-D_DNA_ONLY $(CFLAGS)"


