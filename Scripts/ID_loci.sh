#!/bin/bash -l 

F1=$1
F2=$2
E=$((${F2}*8))


novoindex ${F1}.index ${F1}.fasta 

x=1
while [ $x -le 11 ]
do

	string="sed -n ${x}p data"
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1}')   
	set -- $var
	c1=$1	### Sample name  ###
	
		echo "#!/bin/bash" > ${c1}.sh
		echo "" >> ${c1}.sh

		echo "novoalign  -r E $E -t 180 -d ${F1}.index -f ${F1}.fasta_${c1} > ${F1}.novo_${c1}" >> ${c1}.sh
		
		sbatch -t 1440 --mem=16G -c 1 ${c1}.sh

	x=$(( $x + 1 ))

done


#	cat *novo* > ${F1}.novo
		
#	IdentifyLoci3.pl ${F1}.novo > ${F1}.loci 
#	SimplifyLoci2.pl ${F1}.loci > ${F1}_simple.loci

