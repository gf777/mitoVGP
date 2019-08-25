#!/bin/bash

set -e

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF
	
	This script (map10x2.sh) is used for the final step of short read polishing. This step
	requires the output of the script map10x1.sh to be trimmed at the overlapping ends. This
	can be achieved using mummer or BLAST to determine the coordinates for trimming.
	Experimental script trimmer.sh may also be employed.
	Ends should be trimmed leaving about 100 bp of overlapping ends, in order to achieve a good
	alignment, and those overlapping ends should then be removed from the final assembly.

	It requires the following software (and their dependencies) installed:
	bowtie2/2.3.5, samtools/1.7, freebayes/1.1.0-46-g8d2b3a0-dirty, bcftools/1.9

	Reads are aligned to the reference, Similarly to script map10x1.sh, a final round of
	freebayes and bcftools consensus is required to obtain the polished contig using the 
	aligned outcome of the script (this step is currently not included in the script).

	Required arguments are:
		-s the species name (e.g. Calypte_anna)
		-i the VGP species ID (e.g. bCalAnn1)
		-n the contig ID identified from the BLAST search by the script blastMT.sh
		-t the number of threads

	Optional arguments are:	
		-c if run on cluster. Supported options are:
			SLURM
			None (Default)

EOF

exit 0

fi

printf "\n\n++++ running: map10x2.sh ++++\n\n"

#set options

while getopts ":s:i:n:c:t:" opt; do

	case $opt in
		s)
			SPECIES=$OPTARG
			echo "Species: -s $OPTARG"
			;;
        i)
        	ID=$OPTARG
        	echo "Species ID: -i $OPTARG"
            ;;
		n)
            CONTIG=$OPTARG
			echo "Contig number: -n $OPTARG"
			;;
		c)
            GRID=$OPTARG
			echo "Cluster: -c $OPTARG"
			;;
		t)
			NPROC=$OPTARG
			echo "Number of threads: -t $OPTARG" >&2
            ;;
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac

printf "\n"

done

if [[  ${GRID} == "SLURM" ]]; then

echo Starting at `date`
echo This is job $SLURM_JOB_ID
echo Running on `hostname`

fi

#define working directory
W_URL=${SPECIES}/assembly_MT_rockefeller/intermediates
printf "Working directory: $W_URL\n\n"

if [[ -e "${W_URL}/freebayes_round2/${FNAME}_10x2.fasta" ]]; then

	printf "\n\noutput already present: skipping.\n\n"
	exit 0

fi

if ! [[ -e "${W_URL}/bowtie2_round2" ]]; then

	mkdir ${W_URL}/bowtie2_round2

	FNAME="${ID}.${CONTIG}_arrow2_10x1_trim1"

	#align
	bowtie2-build ${W_URL}/trimmed/${FNAME}.fasta ${W_URL}/bowtie2_round2/${ID}
	bowtie2 -x ${W_URL}/bowtie2_round2/${ID} -1 ${W_URL}/bowtie2_round1/fq/aligned_${ID}_all_1.fq -2 ${W_URL}/bowtie2_round1/fq/aligned_${ID}_all_2.fq -p ${NPROC} --no-mixed | samtools view -bSF4 - > "${W_URL}/bowtie2_round2/aligned_${ID}_all_trimmed.bam"

	#sort and index the alignment
	samtools sort -@ ${NPROC} ${W_URL}/bowtie2_round2/aligned_${ID}_all_trimmed.bam -o ${W_URL}/bowtie2_round2/aligned_${ID}_all_trimmed_sorted.bam -@ ${NPROC}
	samtools index -@ ${NPROC} ${W_URL}/bowtie2_round2/aligned_${ID}_all_trimmed_sorted.bam
	rm ${W_URL}/bowtie2_round2/aligned_${ID}_all_trimmed.bam

fi

if ! [[ -e "${W_URL}/freebayes_round2/" ]]; then

	mkdir ${W_URL}/freebayes_round2/

	freebayes --bam ${W_URL}/bowtie2_round2/aligned_${ID}_all_trimmed_sorted.bam --fasta-reference ${W_URL}/trimmed/${FNAME}.fasta --vcf ${W_URL}/freebayes_round2/aligned_${ID}_all_trimmed_sorted.vcf --theta 0.001 --ploidy 1 --region $(sed -n 1p ${W_URL}/trimmed/${FNAME}.fasta | tr -d '>'):50-$(( $(sed -n 2p ${W_URL}/trimmed/${FNAME}.fasta | tr -d '\n' | wc -c) - 50 ))

	bgzip ${W_URL}/freebayes_round2/aligned_${ID}_all_trimmed_sorted.vcf

	tabix -p vcf ${W_URL}/freebayes_round2/aligned_${ID}_all_trimmed_sorted.vcf.gz

	bcftools consensus ${W_URL}/freebayes_round2/aligned_${ID}_all_trimmed_sorted.vcf.gz -f ${W_URL}/trimmed/${FNAME}.fasta -o ${W_URL}/freebayes_round2/${FNAME}_10x2.fasta

fi