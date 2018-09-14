#!/bin/bash

#Load modules
module load gcccore/6.4.0 bioperl/1.7.1

#Setting directories and files
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
SCRIPTDIR="${WORKDIR}/scripts"
GENEDIR="${WORKDIR}/genes"
GENENAMES=`ls ${GENEDIR}`

while read -r gene; do
echo "processing alignment of gene ${gene}"
ALIGNDIR="${GENEDIR}/${gene}/align"

echo "removing STOP codons"
perl ReplaceStopWithRefCodonGaps.pl  -pep ${ALIGNDIR}/${gene}_align_AA.fasta -nuc ${ALIGNDIR}/${gene}_align_DNA.fasta -output ${ALIGNDIR}/${gene}_align_DNA_non_stop_codons.fasta 

echo "curating alignment"
awk '/^>/{f=!d[$1];d[$1]=1}f' ${ALIGNDIR}/${gene}_align_DNA_non_stop_codons.fasta | sed 's/!/N/g' | cut -d '|' -f 1 | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | awk '{print substr($0, 1, length($0)-3)}' | tr "\t" "\n" > ${ALIGNDIR}/${gene}_align_DNA_curated.fasta

echo "trimming alignment for ${gene}"
trimal -in ${ALIGNDIR}/${gene}_align_DNA_curated.fasta -out ${ALIGNDIR}/${gene}_align_DNA_trimmed2.fasta -gt 0.4 -st 0.1 -resoverlap 0.8 -seqoverlap 70

done <<< "${GENENAMES}"
