Getting Started
---------------

1. Acquire the fermi source code from the [download page][5] and compile with:

       tar -jxf fermi-x.y.tar.bz2
       (cd fermi-x.y; make)

   where `x.y` is the version number.

2. Download the *C. elegans* reads [SRR065390][8] (ftp directory) and convert
   to the FASTQ format with the `fastq-dump` tool from the [SRA toolkit][9]:

       fastq-dump --split-spot SRR065390.lite.sra

3. Perform assembly with:

       fermi-x.y/run-fermi.pl -ct8 -e fermi-x.y/fermi SRR065390.fastq > fmdef.mak
       make -f fmdef.mak -j 8 > fmdef.log 2>&1

   The entire procedure takes about several hours with 8 CPU cores. File
   `fmdef.p5.fq.gz` contains the final contigs. The quality line in the
   FASTQ-like format gives the per-base read depth computed from
   non-redundant error-corrected reads.


FAQ
---

####0. In addition to this FAQ, are there any other documentations?

The algorithms and evaluations are described in the [preprint][1] available
from arXiv. The detailed usage is documented in the [fermi manpage][2].

####1. What is fermi?

Fermi is a de novo assembler for Illumina reads from whole-genome short-gun
sequencing. It also provides tools for error correction, sequence-to-read
alignment and comparison between read sets. It uses the FMD-index, a novel
compressed data structure, as the key data representation.

####2. How is fermi different from other assemblers?

For small genomes, fermi is not much different from other assemblers in terms
of performance. Nonetheless, for mammalian genomes, fermi is one of the few
choices that can do the job in a relatively small memory footprint. It can
assemble 35-fold human data in 90GB shared memory with an overall similar
contiguity and accuracy to other mainstream assemblers.

In addition to de novo assembly, fermi ultimately aims to preserve all the
information in the raw reads, in particular heterozygous events. SNP and INDEL
calling can be achieved by aligning the fermi unitigs to the reference genome
and has been shown to be advantageous over other approaches in some aspects (see
also the [preprint][1]).

####3. What is the relationship between fermi and SGA?

Fermi is substantially influenced by [SGA][3]. It follows a similar workflow,
including the idea of contrasting read sets.  On the other hand, the internal
implementation of fermi is distinct from that of SGA. Fermi is based on a novel
data structure and uses different algorithms for almost every step. As to the
end results, fermi has a similar performance to SGA for features shared between
them, and is arguably easier to use. In all, both fermi and SGA are viable
options for de novo assembly and contrast variant calling.

####4. Are there release notes?

Yes, below this FAQ.

####5. How to install fermi?

You may clone the [fermi github repository][4] to get the latest source code,
or acquire the source code of stable releases from the [download page][5]. You
can compile fermi by invoking `make` in the source code directory. The only
library dependency is [zlib][6]. After compilation, you may copy `fermi` and
`run-fermi.pl` to your `PATH` or simply use the executables in the source code
directory.

####6. How to run fermi for de novo assembly?

The [fermi manpage][2] shows an example. Briefly, if you have Illumina
short-insert paired-end reads `read1.fq.gz` and `read2.fq.gz`, you can run:

    run-fermi.pl -Pe ./fermi -t12 read1.fq.gz read2.fq.gz > fmdef.mak
    make -f fmdef.mak -j 12

to perform assembly using 12 CPU cores. The `fmdef.p5.fq.gz` gives the final
contigs using the paired-end information. If you only want to correct errors,
you may use

    make -f fmdef.mak -j 12 fmdef.ec.fq.gz


[1]: http://arxiv.org/abs/1203.6364
[2]: https://github.com/lh3/fermi/blob/master/fermi.1
[3]: https://github.com/jts/sga
[4]: https://github.com/lh3/fermi
[5]: https://github.com/lh3/fermi/downloads
[6]: http://zlib.net/
[7]: https://github.com/lh3/fermi/blob/master/README.md
[8]: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/litesra/SRX/SRX026/SRX026594/SRR065390/
[9]: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software


Release Notes
-------------

###Release 1.0 (2012-04-09)

This is the first public release of fermi, a de novo assembler and analysis
tool for whole-genome shot-gun sequencing. Source code can be acquired from
the [download page][5]. Please read the [manpage][2] and the [FAQ][7] for
detailed usage.

(1.0: 2012-04-09, r697)
