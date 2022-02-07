#!/bin/bash

#SBATCH --job-name=alignm-map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --time=720

F1=$1

novoindex ${F1}.index ${F1}.fasta 

novoalign  -r E 50 -t 180 -d ${F1}.index -f ${F1}.fasta > ${F1}.novo
		

