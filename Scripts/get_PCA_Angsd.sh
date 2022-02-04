#!/bin/bash

#SBATCH --job-name=pcangsd
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --partition=mid
#SBATCH --mem=120G
#SBATCH --time=20:00:00

module load anaconda/2.7
mkdir results_pca
scripts=~/zonguldak/pop_gen/scripts

##############################################

pop=$1
nInd=$(wc -l ${pop}.bamlist | awk '{print $1}')
mInd=$((${nInd}/2))

#############################################


### Build Geno file in beagle format ###

#	angsd -bam ${pop}.bamlist -out results_pca/${pop} -GL 1 -doMajorMinor 1 -doMaf 1 -doGlf 2 -doGeno 5 -dovcf 1 -doPlink 2 -doPost 1 -postCutoff 0.85 -minMapQ 10 -minQ 20 -minInd $mInd -SNP_pval 1e-12 -minMaf 0.05 -nThreads 2

### Calculate covariance matrix ###

#	python ~/bin/pcangsd/pcangsd.py -beagle results_pca/${pop}.beagle.gz -o results_pca/${pop} -threads 2

### PCA analysis and plot

	Rscript $scripts/plotPCA.R -i results_pca/${pop}.cov -c 1-2 -a ${pop}.clst -o results_pca/${pop}_pca_1_2.pdf
