#!/bin/bash

while IFS=$'\t' read -r species id

do

sbatch --partition=hpc --cpus-per-task=8 mito_qv.sh ${species} ${id} VGP_dataset/${id}*.fasta

done<$1
