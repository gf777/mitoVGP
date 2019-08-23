#!/bin/bash

set -e

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF

	Usage: '$0 -s species -i species_ID -r reference -t threads -l filelist (optional) -c cluster (optional)'

	Required arguments are:
	-s the species name (e.g. Calypte_anna)
	-i the VGP species ID (e.g. bCalAnn1)
	-r the reference sequence fasta file
	-t the number of threads
	
	Optional arguments are:
	-l use files from list of files (default looks into aws)
	-c if run on cluster. Supported options are:
		SLURM
		None (Default)

EOF

exit 0

fi

#set options

printf "\n"

while getopts ":l:s:i:r:g:i:c:o:m:t:w:" opt; do

	case $opt in
		l)
			LIST=$OPTARG
			echo "Mode: -l $OPTARG"
			;;
		s)
			SPECIES=$OPTARG
			echo "Species: -s $OPTARG"
			;;
        i)
        	ID=$OPTARG
        	echo "Species ID: -i $OPTARG"
            ;;
        r)
			REF=$OPTARG
			echo "Reference: -r $OPTARG"
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
W_URL=${SPECIES}/${ID}/assembly_MT_10x

if ! [[ -e "${W_URL}" ]]; then

mkdir -p ${W_URL}

fi

if ! [[ -e "${W_URL}/log" ]]; then

mkdir ${W_URL}/log

dw_date=`date "+%Y%m%d-%H%M%S"`;

fi

if [[ -z ${LIST} ]]; then

	#record 10x raw data files available in the cloud at the time of the analysis
	aws s3 --no-sign-request ls s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/ > ${W_URL}/log/file_list_$dw_date.txt

else

	cat ${LIST} > ${W_URL}/log/file_list_$dw_date.txt

fi

#copy the user-provided reference mitogenome to the reference folder
if ! [[ -e "${W_URL}/reference" ]]; then

	mkdir ${W_URL}/reference
	cp ${REF} ${W_URL}/reference/${REF%.*}.fasta

fi

if ! [[ -e "${W_URL}/${ID}.1.bt2" ]]; then

bowtie2-build ${W_URL}/reference/${REF%.*}.fasta ${W_URL}/reference/${ID}

fi

if [[ -z ${LIST} ]]; then

	#determine fw and rv reads
	if ! [[ -e "${W_URL}/p_fw.txt" ]]; then

		grep -o -E "${ID}.*_R1_.*" ${W_URL}/log/file_list_${dw_date}.txt | sort | uniq > ${W_URL}/p_fw.txt

	fi

	if ! [[ -e "${W_URL}/p_rv.txt" ]]; then

		grep -o -E "${ID}.*_R2_.*" ${W_URL}/log/file_list_${dw_date}.txt | sort | uniq > ${W_URL}/p_rv.txt

	fi
	
else

	grep -o -E ".*_R1_.*" ${W_URL}/log/file_list_${dw_date}.txt | sort | uniq > ${W_URL}/p_fw.txt
	grep -o -E ".*_R2_.*" ${W_URL}/log/file_list_${dw_date}.txt | sort | uniq > ${W_URL}/p_rv.txt

fi

mapfile -t p1 < ${W_URL}/p_fw.txt
mapfile -t p2 < ${W_URL}/p_rv.txt

echo "--Following PE files found:"

for ((i=0; i<${#p1[*]}; i++));
do

echo ${p1[i]} ${p2[i]} $i
printf "\n"

done

if ! [[ -e "${W_URL}/aligned_raw_reads" ]]; then

	mkdir ${W_URL}/aligned_raw_reads

fi

#for each 10x PE raw data do
for ((i=0; i<${#p1[*]}; i++));
do

if ! [[ -e "${W_URL}/aligned_raw_reads/aligned_${ID}_${i}.bam" ]]; then

	echo "--Retrieve and align:"

	if [[ -z ${LIST} ]]; then
	
		echo "s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/${p1[i]}"
		echo "s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/${p2[i]}"
		printf "\n"

		if ! [[ -e "${W_URL}/aligned_raw_reads/${p1[i]}" ]]; then  

		#download R1
		aws s3 --no-sign-request cp s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/${p1[i]} ${W_URL}

		fi

		if ! [[ -e "${W_URL}/${p2[i]}" ]]; then

		#download R2
		aws s3 --no-sign-request cp s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/${p2[i]} ${W_URL}

		fi
	
	else

		if ! [[ -e ${W_URL}/$(basename -- "${p1[i]}") ]]; then
			
			echo "ln -s ${p1[i]} ${W_URL}/$(basename -- "${p1[i]}")"
			ln -s ${p1[i]} ${W_URL}/$(basename -- "${p1[i]}")
		
		fi
			
		if ! [[ -e ${W_URL}/$(basename -- "${p2[i]}") ]]; then
					
			ln -s ${p2[i]} ${W_URL}/$(basename -- "${p2[i]}")
			echo "ln -s ${p2[i]} ${W_URL}/$(basename -- "${p2[i]}")"
			
		fi

		printf "\n"
		
	fi

		p_1=$(basename -- "${p1[i]}")
		p_2=$(basename -- "${p2[i]}")
		
		echo ${W_URL}/${p_1}
		echo ${W_URL}/${p_2}

		if ! [[ -e ${W_URL}/trimmed_${i}_R1_001.fastq.gz ]]; then

		/usr/bin/python bin/proc10xG/process_10xReads.py -1 ${W_URL}/${p_1} -2 ${W_URL}/${p_2} -o ${W_URL}/trimmed_${i}
		
		fi
		
		if [[ -z ${LIST} ]]; then

			rm ${W_URL}/${p1[i]} ${W_URL}/${p2[i]}
		
		fi
		
		#align
		bowtie2 -x ${W_URL}/reference/${ID} -1 ${W_URL}/trimmed_${i}_R1_001.fastq.gz -2 ${W_URL}/trimmed_${i}_R2_001.fastq.gz -p ${NPROC} | samtools view -bSF4 - > "${W_URL}/aligned_raw_reads/aligned_${ID}_${i}.bam"

		#remove
		rm ${W_URL}/trimmed_${i}_R1_001.fastq.gz ${W_URL}/trimmed_${i}_R2_001.fastq.gz

	fi

done

#generate single alignment file out of all raw data
if ! [[ -e "${W_URL}/aligned_${ID}_all.bam" ]]; then

samtools merge ${W_URL}/aligned_${ID}_all.bam ${W_URL}/aligned_raw_reads/*.bam

#downsample the file if >1.5M reads
READ_N=$(bedtools bamtobed -i ${W_URL}/aligned_${ID}_all.bam | cut -f 4 | wc -l)

if (("$READ_N" >= "1500000")); then

N=$(awk "BEGIN {printf \"%.2f\",1500000/${READ_N}}")
N_READ_N=$(awk "BEGIN {printf \"%.0f\",$N*${READ_N}}")

echo "The number of reads is above 1.5M ($READ_N). Downsampling by $N ($N_READ_N)"

samtools view -s $N -b ${W_URL}/aligned_${ID}_all.bam > ${W_URL}/aligned_${ID}_all_sub.bam
mv ${W_URL}/aligned_${ID}_all.bam ${W_URL}/aligned_${ID}_all_o.bam
mv ${W_URL}/aligned_${ID}_all_sub.bam ${W_URL}/aligned_${ID}_all.bam

fi

fi

if ! [[ -e "${W_URL}/fq" ]]; then

mkdir ${W_URL}/fq
samtools fastq ${W_URL}/aligned_${ID}_all.bam -1 ${W_URL}/fq/aligned_${ID}_all_1.fq -2 ${W_URL}/fq/aligned_${ID}_all_2.fq -s ${W_URL}/fq/aligned_${ID}_all_s.fq

fi
