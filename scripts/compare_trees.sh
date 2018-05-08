#!/bin/bash

cd $LUSTRE/cancer_genes_selection
while read -r gene; do
echo /home/uvi/be/avs/lustre/cancer_genes_selection/genes/${gene}/tree/${gene}.raxml.bestTree; done < data/genelist > treelist

module load miniconda
source activate /home/uvi/be/avs/tools

cat treelist | ete3 compare -r species_tree.nwk --src_attr_parser '(^.{6})' --unrooted > tree_distances.txt
