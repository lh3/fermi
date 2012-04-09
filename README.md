####0. In addition to this FAQ, are there any other documentations?

The algorithms and evaluations are described in the
[preprint](http://arxiv.org/abs/1203.6364) available from arXiv. The detailed
command lines are documented in the [fermi
manpage](https://github.com/lh3/fermi/blob/master/fermi.1).

####1. What is fermi?

Fermi is a de novo assembler for Illumina reads from whole-genome short-gun
sequencing. It also provides tools for error correction, sequence-to-read
alignment and comparison between read sets. It uses the FMD-index, a novel
compressed data structure, as the key data representation.

####2. How is fermi different from other assemblers?

For small genomes, fermi is not much different other assemblers in terms of
performance. Nonetheless, for mammalian genomes, fermi is one of the few
choices that can do the job in a relatively small memory footprint. It can
assemble 35-fold human data in 90GB shared memory. It achieves overall similar
contiguity and accuracy to other mainstream assemblers.

In addition to de novo assembly, fermi ultimately aims to preserve all the
information in the raw reads, in particular heterozygous events. SNP and INDEL
calling can be achieved by aligning the fermi unitigs to the reference genome
and has been shown to be advantageous over other approaches to some extend (see
also the [preprint](http://arxiv.org/abs/1203.6364)).

####3. What is the relationship between fermi and SGA?

Fermi is substantially influenced by [SGA](https://github.com/jts/sga). It
follows a similar workflow to SGA, including the idea of contrasting read sets.
On the other hand, the internal implementation of fermi is distinct from that
of SGA. Fermi uses a novel data structure and different algorithms for almost
every step. As to the end results, fermi has a similar performance to SGA in
the intersection of both feature sets and is arguably easier to use. In all,
both fermi and SGA are viable options for de novo assembly and contrast variant
calling.

####4. How to run fermi?

The [fermi manpage](https://github.com/lh3/fermi/blob/master/fermi.1) gives an
example. Briefly, if you have Illumina short-insert paired-end reads `read1.fq.gz`
and `read2.fq.gz`, you can run:

    run-fermi.pl -Pe ./fermi -t12 read1.fq.gz read2.fq.gz > fmdef.mak
    make -f fmdef.mak -j 12

to perform assembly. The manpage explains more details.

