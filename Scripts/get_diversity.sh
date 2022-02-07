#!/bin/bash

#SBATCH --job-name=diversity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=30

##############################################

infile=$1
mkdir results_diversity

#############################################

for i in `cat $infile`;
do

	echo "#!/bin/bash" > ${i}.sh
	echo "" >> ${i}.sh
	echo "realSFS saf2theta results_sfs/${i}.saf.idx -sfs results_sfs/${i}.sfs -outname results_diversity/${i}" >> ${i}.sh
 	echo "thetaStat do_stat results_diversity/${i}.thetas.idx" >> ${i}.sh

	sbatch -J diversity -n 1 -N 1 -p mid -t 720 --mem=32G ${i}.sh

done
