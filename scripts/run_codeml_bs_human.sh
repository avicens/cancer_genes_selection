#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --cpus-per-task 2
#SBATCH -t 05:00:00
#SBATCH --mem 6G

#Setting directories and files
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
GENEDIR="${WORKDIR}/genes"

gene=`ls ${GENEDIR} | sed "${SLURM_ARRAY_TASK_ID}q;d"`

HUMID=`grep -o 'ENSG[0-9]*' ${GENEDIR}/${gene}/tree/${gene}_trimmed.raxml.bestTree` > tmp
echo "the human id for $gene is $HUMID"

module load miniconda
source activate /home/uvi/be/avs/lustre/tools/miniconda_avs

echo "Running Codeml for ${gene}"
ete3 evol -t ${GENEDIR}/${gene}/tree/${gene}_trimmed.raxml.bestTree --alg ${GENEDIR}/${gene}/align/${gene}_align_DNA_trimmed2.fasta --codeml_param CodonFreq,3 ncatG,4 verbose,0 --models bsA bsA1 
 --mark $HUMID --cpu 2 -o ${GENEDIR}/${gene}/paml > ${GENEDIR}/${gene}/paml/${gene}_bs_human.out
