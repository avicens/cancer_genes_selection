#!/bin/env bash

WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
DATADIR="${WORKDIR}/data"
GENEDIR="${WORKDIR}/genes"

while read gene; do
mkdir ${GENEDIR}/$gene
done < ${DATADIR}/genenames.txt
