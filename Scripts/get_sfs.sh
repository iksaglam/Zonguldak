#!/bin/bash

#SBATCH --job-name=run_sfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=30

##############################################

infile=$1
ref=$2
scripts=~/zonguldak/pop_gen/scripts
mkdir results_sfs

#############################################

for i in `cat $infile`;
do

	echo "#!/bin/bash" > ${i}.sh
	echo "" >> ${i}.sh

### Build SAF files ###

	echo "angsd -bam ${i}.bamlist -ref $ref -anc $ref -out results_sfs/$i -GL 1 -doSaf 1 -doCounts 1 -minMapQ 10 -minQ 20" >> ${i}.sh

### Calculate the and plot SFS ###

	echo "realSFS results_sfs/${i}.saf.idx -maxIter 100 > results_sfs/${i}.sfs" >> ${i}.sh
	echo "Rscript $scripts/plotSFS.R results_sfs/${i}.sfs $i 0 results_sfs/${i}.sfs.pdf" >> ${i}.sh

	sbatch -J sfs -n 1 -N 1 -p mid -t 720 --mem=32G ${i}.sh

done
