#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 8
#SBATCH -t 04:00:00
#SBATCH --mem 6G

#Setting directories and files
TOOLSDIR="/home/uvi/be/avs/tools"
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
GENEDIR="${WORKDIR}/genes"

gene=`ls ${GENEDIR} | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Creating directory for tree of ${GENE}"
mkdir ${GENEDIR}/${gene}/tree

module load gcc/5.3.0 openmpi/1.10.2 raxml-ng/0.5.1b

echo "Building tree for ${gene}"
raxml-ng-mpi --all --msa ${GENEDIR}/${gene}/align/${gene}_align_DNA_trimmed.fasta --model GTR+G -tree pars{10} --bs-trees 100 --prefix ${GENEDIR}/${gene}/tree/${gene} --redo
