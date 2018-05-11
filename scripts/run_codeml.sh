#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --cpus-per-task 7
#SBATCH -t 10:00:00
#SBATCH --mem 6G

#Setting directories and files
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
GENEDIR="${WORKDIR}/genes"

gene=`ls ${GENEDIR} | sed "${SLURM_ARRAY_TASK_ID}q;d"`
mkdir ${GENEDIR}/${gene}/paml

module load miniconda
source activate /home/uvi/be/avs/tools/miniconda_avs

echo "Running Codeml for ${gene}"
ete3 evol -t ${GENEDIR}/${gene}/tree/${gene}_trimmed.raxml.bestTree --alg ${GENEDIR}/${gene}/align/${gene}_align_DNA_trimmed2.fasta 
--codeml_param CodonFreq,3 ncatG,4 verbose,0  --models M0 M3 M1 M2 M7 M8 M8a --cpu 7 -o ${GENEDIR}/${gene}/paml > ${GENEDIR}/${gene}/paml/${gene}_test_M7-M8.out
