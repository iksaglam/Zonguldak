#!/bin/bash -l

#SBATCH --job-name=align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=mid
#SBATCH --mem=120G
#SBATCH --time=720

pop=$1
ref=$2
n=$(wc -l ${pop} | awk '{print $1}')
reads=~/zonguldak/pop_gen/data

# bwa index -a bwtsw ${ref}
	
for i in `cat $pop`;
do

		echo "#!/bin/bash" > ${i}.sh
		echo "" >> ${i}.sh
		echo "bwa mem ${ref} $reads/${i}_R1.fastq.gz $reads/${i}_R2.fastq.gz > ${i}.sam" >> ${i}.sh
		echo "samtools view -bS ${i}.sam > ${i}.bam" >> ${i}.sh
		echo "samtools sort ${i}.bam ${i}_sorted" >> ${i}.sh
		echo "samtools view -b -f 0x2 ${i}_sorted.bam > ${i}_sorted_proper.bam" >> ${i}.sh
		echo "java -jar /kuacc/apps/picard/2.22.1/picard.jar MarkDuplicates INPUT=${i}_sorted_proper.bam OUTPUT=${i}_sorted_proper_rmdup.bam METRICS_FILE=${i}_metrics.txt VALIDATION_STRINGENCY=LENIENT  REMOVE_DUPLICATES=True" >> ${i}.sh
		echo "samtools index ${i}_sorted_proper_rmdup.bam" >> ${i}.sh

		sbatch -J iks-align -n 1 -N 1 -p mid -t 720 --mem=32G ${i}.sh

done

