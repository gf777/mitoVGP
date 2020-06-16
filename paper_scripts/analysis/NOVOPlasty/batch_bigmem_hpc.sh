#!/bin/bash

while IFS=$'\t' read -r species id ref

do

sbatch --partition=bigmem --cpus-per-task=21 NOVOplasty.sh $species $id *${ref}*

done < $1
