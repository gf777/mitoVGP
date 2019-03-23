#!/bin/bash

#this script (mtDNApipe.sh) is used to retrieve mitochondrial-like sequences from the raw Pacbio data 
#generated in the framework of the Vertebrate Genomes Project and assemble them using Canu

#it requires the following dependencies installed:
#blasR, bedtools2 (bamToFastq), Canu, blastn

#sequence retrieval is based on a search by similarity using BLASR alignment
#all Pacbio raw data files are first downloaded from the Genomeark and then individually
#aligned to a reference genome provided by the user
#the reference genome can be from the same species if available, or from a
#closely-to-distantly related species
#the approach is similar to that of Organelle_PBA described in:
#Soorni et al. BMC Genomics (2017) DOI 10.1186/s12864-016-3412-9

#in the second steps reads are used to generate assemblies using Canu assembler
#usually with default parameters
#(for rGopEvg1 and mRhiFer1 minOverlapLength=300 correctedErrorRate=0.105 were used, respectively)

#the reference genome provided by the user is then blasted to the contigs generated 
#by Canu to identify the putative mitocontig

#required positional arguments are:
#the species name (e.g. Calypte_anna)
#the VGP species ID (e.g. bCalAnn1)
#the putative mitogenome size (potentially, that of the reference genome). It does not
#need to be precise

set -e

#set variables species abbreviation
SPECIES=$1
ABBR=$2
SIZE=$3

dw_date=`date "+%Y%m%d-%H%M%S"`; aws s3 --no-sign-request ls s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/pacbio/ | grep -oP "m.*.subreads.bam" | uniq > ${SPECIES}/file_list_$dw_date.txt

while read p; do

if ! [[ $p == *scraps* ]] && ! [[ $p == *.pbi ]] && [[ $p == *.bam ]] && ! [[ -e "${SPECIES}/aligned_${p}" ]] && ! [[ -e "${SPECIES}/aligned_pb_raw_reads/aligned_${p}" ]]; then

aws s3 --no-sign-request cp s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/pacbio/$p ${SPECIES}
blasr --bestn 1 ${SPECIES}/$p ${SPECIES}/mtDNA_${SPECIES}.fasta --bam --out ${SPECIES}/aligned_$p --nproc 16
rm ${SPECIES}/$p

fi

done <${SPECIES}/file_list_$dw_date.txt

cd ${SPECIES}

if ! [[ -e "aligned_pb_raw_reads" ]]; then

mkdir aligned_pb_raw_reads

fi

if [ -f aligned_*.bam ]; then

mv aligned_*.bam aligned_pb_raw_reads/

fi

if ! [[ -e "fq" ]]; then

#convert to fq and merge into a single read file
for f in aligned_pb_raw_reads/aligned_*.bam; do filename=$(basename -- "$f"); filename="${filename%.*}"; bamToFastq -i $f -fq "${filename}.fq"; done

mkdir fq
mv *.fq fq/

fi


if ! [[ -e "${ABBR}.fastq" ]]; then

cat fq/*.fq > ${ABBR}.fastq

fi

if ! [[ -e "${ABBR}_canu" ]]; then

#assemble mtDNA reads with canu
canu \
 -p ${ABBR} -d ${ABBR}_canu \
 genomeSize=${SIZE} \
 -pacbio-raw ${ABBR}.fastq \
 correctedErrorRate=0.105

fi

mkdir -p ${ABBR}_blast
cd ${ABBR}_blast

#build blast db
makeblastdb -in ../${ABBR}_canu/${ABBR}.contigs.fasta -parse_seqids -dbtype nucl -out ${ABBR}.db

#blastn
blastn -query ../mtDNA_${SPECIES}.fasta -db ${ABBR}.db -out ${ABBR}.out
blastn -query ../mtDNA_${SPECIES}.fasta -db ${ABBR}.db -outfmt 6 -out ${ABBR}.tb
sed -i "1iquery_acc.ver\tsubject_acc.ver\t%_identity\talignment_length\tmismatches\tgap_opens\tq.start\tq.end\ts.start\ts.end\tevalue\tbitscore" ${ABBR}.tb
cat ${ABBR}.tb | column -t > ${ABBR}_results.txt