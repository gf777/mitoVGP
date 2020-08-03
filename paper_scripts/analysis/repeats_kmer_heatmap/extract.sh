#!/bin/bash

while IFS=$'\t' read -r c1 c2 c3 c4 c5

do

if [[ ! -e regions/${c2}_${c5}.fasta ]]; then

samtools faidx VGP_dataset/${c2}*
samtools faidx VGP_dataset/${c2}*.fasta ${c2}:${c3}-${c4} > regions/${c2}_${c5}.fasta

fi

if [[ ! -d meryl/${c2}_${c5}_k${2} ]]; then

meryl count k=${2} regions/${c2}_${c5}.fasta output meryl/${c2}_${c5}_k${2}

fi

done<$1

f=($(find meryl -maxdepth 1 -type d -printf '%P\n' | grep "k${2}"))

truncate -s 0 results${2}.txt

for ((i = 0; i < ${#f[@]}; i++)); do 
      for ((j = i + 1; j < ${#f[@]}; j++)); do 

		a=$(meryl print meryl/${f[i]} 2> /dev/null | awk '{sum+=$2}END{print sum}')
		b=$(meryl print meryl/${f[j]} 2> /dev/null | awk '{sum+=$2}END{print sum}')
		c=$(meryl print intersect-sum meryl/${f[i]} meryl/${f[j]} 2> /dev/null | awk '{sum+=$2}END{print sum}')

		if [ -z "$c" ]; then
			c=0
		fi

		d=$(echo "$c / ($a + $b)" | bc -l)

		printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${f[i]}" "${f[j]}" $a $b $c $d
		
		printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${f[i]}" "${f[j]}" $d >> results${2}.txt

	done;
done

awk '{printf $2"\t"$1"\t"$3"\n"}' results${2}.txt > results_rev${2}.txt

cat results${2}.txt results_rev${2}.txt > results_combined${2}.txt

datamash -s crosstab 2,1 unique 3 < results_combined${2}.txt > matrix${2}.txt 
    
sed -i 's/N\/A/1/g' matrix${2}.txt 

rm results${2}.txt results_rev${2}.txt results_combined${2}.txt
