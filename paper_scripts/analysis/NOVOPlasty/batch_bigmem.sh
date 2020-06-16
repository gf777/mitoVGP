#!/bin/bash

while IFS=$'\t' read -r species id ref

do

sbatch --partition=vgl_bigmem --cpus-per-task=32 NOVOplasty.sh $species $id *${ref}*

done < $1
