#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 00:05:00
#SBATCH --mem 6G

#Setting directories and files
TOOLSDIR="/home/uvi/be/avs/lustre/tools"
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
GENEDIR="${WORKDIR}/genes"

GENE=`ls ${GENEDIR} | sed "1q;d"`
echo "Creating directory for alignment of ${GENE}"
mkdir ${GENEDIR}/${GENE}/align
seqfile="${GENEDIR}/${GENE}/seqs/${GENE}_seqs.fas"


module load jdk/1.8.0

echo "Running codon alignment of ${GENE} sequences"
java -jar ${TOOLSDIR}/macse_v0.9b1.jar -i ${seqfile} -g -7 -x -1 -f -30 -d 1 -s -100 -o ${GENEDIR}/${GENE}/align/cds/${GENE}_cds
