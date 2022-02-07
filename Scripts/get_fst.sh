#!/bin/bash

#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=30

##############################################

infile=$1
n=$(wc -l $infile | awk '{print $1}')
mkdir results_fst

#############################################

x=1
while [ $x -le $n ] 
do
	y=$(( $x + 1 ))
	while [ $y -le $n ]
	do

	pop1=$( (sed -n ${x}p $infile) )  
	pop2=$( (sed -n ${y}p $infile) )

	echo "#!/bin/bash" > ${pop1}.${pop2}.sh
	echo "" >> ${pop1}.${pop2}.sh
	echo "realSFS fst index results_sfs/${pop1}.saf.idx results_sfs/${pop2}.saf.idx -sfs results_sfs/${pop1}.${pop2}.2dsfs -fstout results_fst/${pop1}.${pop2}" >> ${pop1}.${pop2}.sh
	echo "realSFS fst stats results_fst/${pop1}.${pop2}.fst.idx 2> results_fst/${pop1}.${pop2}_global.fst" >> ${pop1}.${pop2}.sh
 
	sbatch -J fst -n 1 -N 1 -p mid -t 720 --mem=32G ${pop1}.${pop2}.sh
 
	y=$(( $y + 1 ))   

	done

x=$(( $x + 1 ))

done
