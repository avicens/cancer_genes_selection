#!/bin/bash

#Setting directories and files
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
GENEDIR="${WORKDIR}/genes"
GENENAMES=`ls ${GENEDIR}`

while read -r gene; do
echo "trimming alignment for ${gene}"
trimal -in ${GENEDIR}/${gene}/align/${gene}_align_DNA_trimmed.fasta -out ${GENEDIR}/${gene}/align/${gene}_align_DNA_trimmed2.fasta -gt 0.4 -st 0.1 -resoverlap 0.8 -seqoverlap 70
done <<< "${GENENAMES}"
