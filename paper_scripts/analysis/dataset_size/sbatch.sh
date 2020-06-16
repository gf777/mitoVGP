#!/bin/bash

while IFS=$'\t' read -r species id
do

sbatch --partition=hpc,vgl,vgl_bigmem --cpus-per-task=8 download.sh ${species} ${id}

done<$1
