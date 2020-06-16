#!/bin/bash

while read bam; do

sbatch --partition=vgl,hpc bam2fastq.sh ${bam}

done<${1}
