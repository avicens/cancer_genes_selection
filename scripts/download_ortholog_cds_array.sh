#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -t 00:30:00
#SBATCH --mem 6G

#Setting directories and files
WORKDIR="/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection"
DATADIR="${WORKDIR}/data"
GENEDIR="${WORKDIR}/genes"
DBDIR="${WORKDIR}/db/pep"
SCRIPTDIR="${WORKDIR}/scripts"
SPNAMES=${DATADIR}/species_codes.txt

GENE=`ls ${GENEDIR} | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Creating directory for coding sequences of ${GENE}"
mkdir ${GENEDIR}/${GENE}/seqs
protlist="${GENEDIR}/${GENE}/*_orthologs2.txt"

#Load PERL module
module load perl/5.24.0
#Set path to Biomart modules
export PERL5LIB=$PERL5LIB:/home/uvi/be/avs/tools/biomart-perl/lib/:/home/uvi/be/avs/tools/biomart-perl/perl_modules/

while read -r id; do
spcode=`echo $id | cut -c -6`
dataset=`grep $spcode ${SPNAMES} | cut -d' ' -f2`
echo "Downloading ortholog coding sequences of ${GENE} from ${dataset}"
perl ${SCRIPTDIR}/download_cds_from_protid.pl $id $dataset >> ${GENEDIR}/${GENE}/seqs/"${GENE}_seqs.fas"
done < ${protlist}
