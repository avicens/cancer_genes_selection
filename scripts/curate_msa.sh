#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 20:00:00
#SBATCH --mem 6G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=avicens@uvigo.es

#Setting directories and files
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
GENEDIR="${WORKDIR}/genes"
GENENAMES=`ls ${GENEDIR}`

while read -r gene; do
echo "curating alignment for ${gene}"
cd ${GENEDIR}/${gene}/align/
awk '/^>/{f=!d[$1];d[$1]=1}f' ${gene}_align_DNA.fasta | sed 's/!/N/g' | cut -d'|' -f1 > ${gene}_align_DNA_curated.fasta

done <<< "${GENENAMES}"
