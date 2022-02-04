#!/bin/bash -l

#SBATCH --job-name=admix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=30

##############################################

infile=$1 ### pop name ##
k=$2  ### k number ###
mkdir results_admix

#############################################


### Calculate admixture proportions and likelihoods ###
 
x=2
while [ $x -le $k ] 
do
	y=1
	while [ $y -le 15 ]
	do

	echo "#!/bin/bash" > admix${x}_run${y}.sh
	echo "" >> admix${x}_run${y}.sh
	echo "~/bin/NGSadmix -likes results_pca/${infile}.beagle.gz -K $x -P 4 -seed $[RANDOM] -o results_admix/${infile}_admix${x}_run${y}" >> admix${x}_run${y}.sh

	sbatch -J admix -p mid  -t 720 -n 1 -N 1 admix${x}_run${y}.sh

	y=$(( $y + 1 ))

	done

x=$(( $x + 1 ))
 
done

##(for log in `ls *.log`; do grep -Po 'like=\K[^ ]+' $log; done) > logfile
##(for log in `ls *.log`; do grep -Po 'nPop=\K[^ ]' $log; done) > logfile

