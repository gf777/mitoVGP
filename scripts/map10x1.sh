#!/bin/bash

set -e

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF

	This script (map10x1.sh) is used for short read polishing of the long-read assembly 
	resulting from Canu (mtDNApipe.sh) and Arrow polishing (ArrowPolish.sh).
	In the VGP pipeline it uses 10x data (currently with random barcodes) but it can be used
	with any short read data containing mitogenomic reads using information from Canu
	intermediate files.

	It requires the following software (and their dependencies) installed:
	bowtie2/2.3.5, aws-cli/1.16.101, samtools/1.7, freebayes/1.1.0-46-g8d2b3a0-dirty, bcftools/1.9

	Reads are aligned to the reference. A final round of freebayes and bcftools consensus
	is required to obtain the polished contig using the aligned outcome of the script (this step is
	currently not included in the script).

	In addition, the script provides the fw and rv reads (extracted from the alignment)
	required for the final polishing step (map10x2.sh).

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

printf "\n\n++++ running: map10x1.sh ++++\n\n"

if [[ -e "${W_URL}/freebayes_round1/${ID}.${CONTIG}_arrow2_10x1.fasta" ]]; then

	printf "\n\noutput already present: skipping.\n\n"
	exit 0

fi

#set options

while getopts ":s:i:n:c:t:d:" opt; do

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
		d)
			DOWNL=$OPTARG
			echo "Multithreaded download: -d $OPTARG"
            ;;			
	esac
	
printf "\n"

done

if [[  ${GRID} == "SLURM" ]]; then

	echo Starting at `date`
	echo This is job $SLURM_JOB_ID
	echo Running on `hostname`
	printf "\n"

fi


#define working directory
W_URL=${SPECIES}/assembly_MT_rockefeller/intermediates
printf "Working directory: $W_URL\n\n"

if ! [[ -e "${W_URL}/bowtie2_round1" ]]; then

	mkdir ${W_URL}/bowtie2_round1

fi

#record 10x raw data files available in the cloud at the time of the analysis

dw_date=`date "+%Y%m%d-%H%M%S"`

if ! [[ -e "${W_URL}/bowtie2_round1/log/" ]]; then

	mkdir ${W_URL}/bowtie2_round1/log

fi

if ! [[ -e "${W_URL}/bowtie2_round1/log/file_list_$dw_date.txt" ]]; then

	aws s3 ls s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/ > ${W_URL}/bowtie2_round1/log/file_list_$dw_date.txt

fi

if ! [[ -e "${W_URL}/bowtie2_round1/${ID}.1.bt2" ]]; then

bowtie2-build ${W_URL}/arrow/arrow_round2/${ID}.${CONTIG}_arrow2.fasta ${W_URL}/bowtie2_round1/${ID}

fi

#determine fw and rv reads
if ! [[ -e "${W_URL}/bowtie2_round1/p_fw.txt" ]]; then

	grep -o -E "${ID}.*_R1_.*" ${W_URL}/bowtie2_round1/log/file_list_${dw_date}.txt | sort | uniq > ${W_URL}/bowtie2_round1/p_fw.txt

fi

if ! [[ -e "${W_URL}/bowtie2_round1/p_rv.txt" ]]; then

	grep -o -E "${ID}.*_R2_.*" ${W_URL}/bowtie2_round1/log/file_list_${dw_date}.txt | sort | uniq > ${W_URL}/bowtie2_round1/p_rv.txt

fi

mapfile -t p1 < ${W_URL}/bowtie2_round1/p_fw.txt
mapfile -t p2 < ${W_URL}/bowtie2_round1/p_rv.txt

echo "\n--Following PE files found:"

for ((i=0; i<${#p1[*]}; i++));
do

echo ${p1[i]} ${p2[i]} $i

done

printf "\n"

if ! [[ -e "${W_URL}/bowtie2_round1/aligned_raw_reads" ]]; then

	mkdir ${W_URL}/bowtie2_round1/aligned_raw_reads

fi

if [[ ${DOWNL} == true ]] && ! [[ "$(ls -A ${W_URL}/bowtie2_round1/aligned_raw_reads)" ]] ; then

	aws s3 cp --recursive --include="*.fastq.gz" --exclude="*I1*" s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/ ${W_URL}

fi

#for each 10x PE raw data do
for ((i=0; i<${#p1[*]}; i++));
do

	if ! [[ -e "${W_URL}/bowtie2_round1/aligned_raw_reads/aligned_${ID}_${i}.bam" ]]; then

		printf "\n--Align:\n"

		printf "${p1[i]}"
		printf "${p2[i]}"

		if [[ -z ${DOWNL} ]] || ! [[  ${DOWNL} == true ]]; then
			#download
			aws s3 cp s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/${p1[i]} ${W_URL}
			aws s3 cp s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/${p2[i]} ${W_URL}
		fi
		
		#align
		bowtie2 -x ${W_URL}/bowtie2_round1/${ID} -1 ${W_URL}/${p1[i]} -2 ${W_URL}/${p2[i]} -p ${NPROC} | samtools view -bSF4 - > "${W_URL}/bowtie2_round1/aligned_raw_reads/aligned_${ID}_${i}.bam"

		#remove
		rm ${W_URL}/${p1[i]} ${W_URL}/${p2[i]}

	fi

done

#generate single alignment file out of all raw data
if ! [[ -e "${W_URL}/bowtie2_round1/aligned_${ID}_all.bam" ]]; then

	samtools merge ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam ${W_URL}/bowtie2_round1/aligned_raw_reads/*.bam

	#downsample the file if >1.5M reads
	#READ_N=$(bedtools bamtobed -i ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam | cut -f 4 | wc -l)

	# if (("$READ_N" >= "1500000")); then
	# 
	# N=$(awk "BEGIN {printf \"%.2f\",1500000/${READ_N}}")
	# N_READ_N=$(awk "BEGIN {printf \"%.0f\",$N*${READ_N}}")
	# 
	# echo "The number of reads is above 1.5M ($READ_N). Downsampling by $N ($N_READ_N)"
	# 
	# samtools view -s $N -b ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam > ${W_URL}/bowtie2_round1/aligned_${ID}_all_sub.bam
	# mv ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam ${W_URL}/bowtie2_round1/aligned_${ID}_all_o.bam
	# mv ${W_URL}/bowtie2_round1/aligned_${ID}_all_sub.bam ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam
	# 
	# fi

fi

#split the fw and rv reads in the alignment for next step (map10x2.sh)
if ! [[ -e "${W_URL}/bowtie2_round1/fq" ]]; then

	mkdir ${W_URL}/bowtie2_round1/fq
	samtools fastq ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam -1 ${W_URL}/bowtie2_round1/fq/aligned_${ID}_all_1.fq -2 ${W_URL}/bowtie2_round1/fq/aligned_${ID}_all_2.fq -s ${W_URL}/bowtie2_round1/fq/aligned_${ID}_all_s.fq

fi

#sort and index the alignment
if ! [[ -e "${W_URL}/bowtie2_round1/aligned_${ID}_all_sorted.bam" ]]; then

	samtools sort ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam -o ${W_URL}/bowtie2_round1/aligned_${ID}_all_sorted.bam -@ ${NPROC}
	samtools index ${W_URL}/bowtie2_round1/aligned_${ID}_all_sorted.bam

rm ${W_URL}/bowtie2_round1/aligned_${ID}_all.bam

fi

if ! [[ -e "${W_URL}/freebayes_round1/" ]]; then

	mkdir ${W_URL}/freebayes_round1/

fi

if ! [[ -e "${W_URL}/freebayes_round1/aligned_${ID}_all_sorted.vcf" ]] && ! [[ -e "${W_URL}/freebayes_round1/aligned_${ID}_all_sorted.vcf.gz" ]]; then

	freebayes -f ${W_URL}/arrow/arrow_round2/${ID}.${CONTIG}_arrow2.fasta -b ${W_URL}/bowtie2_round1/aligned_${ID}_all_sorted.bam --ploidy 1 -v ${W_URL}/freebayes_round1/aligned_${ID}_all_sorted.vcf

fi

if ! [[ -e "${W_URL}/freebayes_round1/aligned_${ID}_all_sorted.vcf.gz" ]]; then

	bgzip ${W_URL}/freebayes_round1/aligned_${ID}_all_sorted.vcf -@ ${NPROC}

fi

if ! [[ -e "${W_URL}/freebayes_round1/${ID}.${CONTIG}_arrow2_10x1.fasta" ]]; then

	tabix -p vcf ${W_URL}/freebayes_round1/aligned_${ID}_all_sorted.vcf.gz

	bcftools consensus ${W_URL}/freebayes_round1/aligned_${ID}_all_sorted.vcf.gz -f ${W_URL}/arrow/arrow_round2/${ID}.${CONTIG}_arrow2.fasta -o ${W_URL}/freebayes_round1/${ID}.${CONTIG}_arrow2_10x1.fasta

fi