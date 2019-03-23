#!/bin/bash

set -e

SPECIES=$1
ABBR=$2
CONTIG=$3

cd ${SPECIES}/${ABBR}_10x2

bowtie2-build ${CONTIG} ${ABBR}
bowtie2 -x ${ABBR} -1 ../${ABBR}_10x1/fq/aligned_${ABBR}_all_1.fq -2 ../${ABBR}_10x1/fq/aligned_${ABBR}_all_2.fq -p 16 --no-mixed | samtools view -bSF4 - > "aligned_${ABBR}_all_trimmed.bam"

samtools sort aligned_${ABBR}_all_trimmed.bam -o aligned_${ABBR}_all_trimmed_sorted.bam -@ 16
samtools index aligned_${ABBR}_all_trimmed_sorted.bam