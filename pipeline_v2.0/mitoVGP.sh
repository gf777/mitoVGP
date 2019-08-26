#!/bin/bash

set -e

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF

	Usage: '$0 -s species -i species_ID -r reference -g genome_size -t threads -m mapper (optional) -l filelist (optional) -f size_cutoff (optional) -o canu_options (optional)'

	mitoVGP is used for de novo mitogenome assembly using a combination of long and short read data.
	
	Che the github page https://github.com/GiulioF1/mitoVGP for a description of the pipeline.
	
	This script a simple wrapper of the scripts found in the scripts/ folder. You can find more information
	on each step in the help (-h) of each script.

	Required arguments are:
	-s the species name (e.g. Calypte_anna)
	-i the VGP species ID (e.g. bCalAnn1)
	-r the reference sequence fasta file
	-g the putative mitogenome size (potentially, that of the reference genome). It does not
	need to be precise. Accepts Canu formatting.
	-t the number of threads
	
	Optional arguments are:
	-l use files from list of files (default looks into aws)
	-m the aligner (blasr|minimap2|pbmm2). Default is pbmm2
	-f filter reads by size prior to assembly (reduces the number of NUMT reads and helps the assembly)
	-o the options for Canu

EOF

exit 0

fi

#set options

printf "\n"

while getopts ":s:i:r:g:t:l:m:f:o" opt; do

	case $opt in
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
		g)
			SIZE=$OPTARG
			echo "Genome size: -g $OPTARG"
            ;;
		t)
			NPROC=$OPTARG
			echo "Number of threads: -t $OPTARG"
            ;;
		l)
			LIST=$OPTARG
			echo "Mode: -l $OPTARG"
			;;
		m)
			ALN=$OPTARG
			echo "Aligner: -m $OPTARG"
			;;
		f)
            FL=$OPTARG
			echo "Canu options: -f $OPTARG"
			;;
		o)
            OPTS=$OPTARG
			echo "Canu options: -o $OPTARG"
			;;
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac

printf "\n"

done

#define working directory
W_URL=${SPECIES}/assembly_MT_rockefeller/intermediates

if ! [[ -e "${W_URL}" ]]; then

	mkdir -p ${W_URL}

fi

if ! [[ -e "${W_URL}/log" ]]; then

	mkdir ${W_URL}/log

fi

#copy the user-provided reference mitogenome to the reference folder
if ! [[ -e "${W_URL}/reference" ]]; then

	mkdir ${W_URL}/reference
	cp ${REF} ${W_URL}/reference/${REF%.*}.fasta

fi

#retrieve mito-like reads and assemble
scripts/mtDNApipe.sh -s ${SPECIES} -i ${ID} -r ${REF} -g ${SIZE} -t ${NPROC} 2>&1 | tee ${W_URL}/log/${ID}_mtDNApipe.out &&

#identify the mitocontig
scripts/blastMT.sh -s ${SPECIES} -i ${ID} -r ${REF} 2>&1 | tee ${W_URL}/log/${ID}_blastMT.out &&

CONTIG_ID=$(cat ${W_URL}/blast/${ID%.*.*}_candidate_mitocontig.txt) &&

#polish the mitocontig with long reads
scripts/ArrowPolish.sh -s ${SPECIES} -i ${ID} -n ${CONTIG_ID} -t ${NPROC} 2>&1 | tee ${W_URL}/log/${ID}_ArrowPolish.out &&

#polish the mitocontig with short reads
scripts/map10x1.sh -s ${SPECIES} -i ${ID} -n ${CONTIG_ID} -t ${NPROC} 2>&1 | tee ${W_URL}/log/${ID}_map10x1.out &&

#trim the mitocontig
scripts/trimmer.sh -s ${SPECIES} -i ${ID} -n ${CONTIG_ID} -t ${NPROC} 2>&1 | tee ${W_URL}/log/${ID}_trimmer.out &&

#polish the trimmed mitocontig with short reads
scripts/map10x2.sh -s ${SPECIES} -i ${ID} -n ${CONTIG_ID} -t ${NPROC} 2>&1 | tee ${W_URL}/log/${ID}_map10x2.out &&

#perform final end trimming
scripts/trimmer2.sh -s ${SPECIES} -i ${ID} -n ${CONTIG_ID} -t ${NPROC} 2>&1 | tee ${W_URL}/log/${ID}_trimmer2.out &&