#!/bin/bash

while IFS=$'\t' read -r species id ref

do

sbatch --nice=10000 --partition=vgl --cpus-per-task=32 NOVOplasty.sh $species $id refs/*${ref}*

done < $1
