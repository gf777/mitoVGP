#!/bin/bash

echo Starting at `date`
echo This is job $SLURM_JOB_ID
echo Running on `hostname`

/rugpfs/fs0/vgl/store/vglshare/tools/VGP-tools/canu-SLURM-Rockefeller/Linux-amd64/bin/canu -genomeSize=17508 -p bCucCan1 -d bCuCan1_filtered -pacbio-raw filtered_bCucCan1.fastq

