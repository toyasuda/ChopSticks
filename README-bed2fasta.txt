This document explains how to use scripts that generates simulated genomes.


Usage: perl bed2fasta.pl [options] <FASTA file> <BED file>


DESCRIPTIONS

   Generates a simulated genome in FASTA format based on given <FASTA file> that
   contains SVs specified with <BED file> and SNPs randomly generated
   at rate <SNP rate>.
   Here are options available:

   -r<seed>
       Specifies the seed of random numbers to be genrated. By giving the
       seed explicitly with this option, you can make bed2fasta.pl generates
       exactly the same simulated genome whenever you execute this script.

   -s<SNP rate>
       Generates randome SNPs at each base with probability given by this option.

   -w<number of bases in each line>
       Generates SNPs at each base with probability given by this option.

EXAMPLE

   $ cat lines2.fa
   >genome
   ACTCGAATTCCGAATAGATAG
   $ cat lines2-sv.bed
   genome  1       1       INS     TTTTTTTTTT
   genome  4       8       DEL     4bp
   $ perl ../../ChopSticks/bed2fasta.pl -r1 -s0.01 lines2.fa lines2-sv.bed
   # Using seed 1.
   # Reading genome
   # Reading done. (last: genome)
   genome	5	6	SNP	A	T
   genome	20	21	SNP	G	T
   genome  1       1       INS     TTTTTTTTTT
   genome  4       8       DEL     4bp
   >genome
   ATTTTTTTTTTCTCTCCGAATAGATAT
   $

In the example above, SNP at base 5 is ignored because of deletion
specified by lines2-sv.bed.


COPYRIGHT
      Copyright (C) 2012 Tomohiro Yasuda.
      License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
      This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

AUTHOR
      Tomohiro Yasuda
      The Human Genome Center, Institute of Medical Science, the University of Tokyo.
      t y a s u d a (at) h g c (dot) j p

