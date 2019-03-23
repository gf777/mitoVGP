#!/bin/bash

set -e

#set variables species abbreviation
SPECIES=$1

for f in *.bam; do

if ! [[ $f == *scraps* ]] && ! [[ $f == *.pbi ]] && [[ $f == *.bam ]]; then

 ~/bin/miniconda3/bin/blasr $f ../mtDNA_$SPECIES.fasta --bam --out aligned_$f --nproc 16
 
fi

done

mkdir aligned_raw_reads
mv aligned_*.bam aligned_raw_reads/

mkdir raw_reads
mv *.bam raw_reads/