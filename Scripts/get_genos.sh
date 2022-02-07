#!/bin/bash

#SBATCH --job-name=genotype
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --partition=mid
#SBATCH --mem=120G
#SBATCH --time=20:00:00


##############################################

pop=$1
ref=$2
nInd=$(wc -l ${pop}.bamlist | awk '{print $1}')
mInd=$((${nInd}/2))
mkdir results_geno

#############################################


### Build Geno file in beagle format ###
	angsd -bam ${pop}.bamlist -ref $ref -out results_geno/${pop} -GL 1 -doMajorMinor 1 -doMaf 1 -doGlf 2 -doGeno 5 -dovcf 1 -doPlink 2 -doPost 1 -postCutoff 0.85 -minMapQ 10 -minQ 20 -minInd $mInd -SNP_pval 1e-12 -minMaf 0.05 -nThreads 2

