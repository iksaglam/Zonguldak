#!/bin/bash

F1=$1
data=/kuacc/users/iksaglam/zonguldak/id_rad_loci/data

for i in `cat $F1`;
do

echo "#!/bin/bash" > ${i}.sh
echo "" >> ${i}.sh
echo "perl ./QualityFilter.pl $data/${i}_R1.fastq > ${i}_L80P80.fastq" >> ${i}.sh
echo "perl ./HashSeqs.pl ${i}_L80P80.fastq ${i} > ${i}_L80P80.hash" >> ${i}.sh

sbatch -t 30 ${i}.sh

done

