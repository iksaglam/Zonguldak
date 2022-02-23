#!/bin/bash -l

#SBATCH --job-name=count
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=mid
#SBATCH --time=720

pop=$1
reads=~/zonguldak/pop_gen/data



for i in `cat ${pop}.list`;
do
	
		samtools view -c -f 3 -F 12 ${i}_sorted.bam >> no_align.txt
		samtools view -c -f 1 -F 12 ${i}_sorted_proper.bam >> no_prop.txt
		samtools view -c -f 1 -F 12 ${i}_sorted_proper_rmdup.bam >> no_rmdup_align.txt
		zless $reads/${i}_R1.fastq.gz | wc -l >> no_R1_reads.txt
			

done

	awk '{c=$1/2; print c}' no_R1_reads.txt > no_reads.txt
	awk '{print $1}' ${pop}.list > names.txt
	paste names.txt no_reads.txt no_align.txt no_prop.txt no_rmdup_align.txt > ${pop}.reads
	sed -i '1iIndv	no_reads	no_align	no_prop	no_rmdup' ${pop}.reads
	rm *.txt


