# Molecular Adaption of Cancer Genes
## Description
In this repository I describe the workflow followed to process and analyze the data of the study "Selective constraints driving the Germinal Evolution of Cancer Genes in Mammals" *(manuscript in preparation)*.
The repository contains the scripts and some example data used in the study.

## Cancer genes dataset
I retrieved cancer genes from the Cancer Gene Census (CGC) repository of the COSMIC v84 database. A list of 574 genes classified as Tier 1 (i.e. those genes with a documented activity in cancer) was downloaded. This dataset can be accessed with the link https://cancer.sanger.ac.uk/cosmic/census?tier=1 on March 5th, 2018. The dataset correspons to the **Table S1** of the manuscript. 

The dataset was cleaned in order to be used for later analysis with the R script **build_fundata_table.R** 

The variables included in the dataset are:
* Somatic
* Germline
* Tumor type
* Tissue type
* Molecular genetics
* Role in cancer
* Mutation type

## Retrieving human protein annotations
The first column of the cancer genes dataset, which contains the gene names, was extracted and the information for human genes was retrieved using Ensembl BioMArt (http://www.ensembl.org/biomart/martview/4ee102879139fc3bf745f9a867064956). In the exported file, genes with no protein annotated were discarded, obtaining 535 genes. On these genes, the protein isoform with the best transcript support level was chosen. The obtained dataset is the file **ensembl_gene_uniq.tsv**.

In this step, and in order to apply the later analysis per gene, a working folder was created for each gene with the bash script **create_gene_folders.sh**.

The Ensembl_protein_ID column was extracted from the *ensembl_gene_uniq.tsv* table, getting the **protein_ensembl_id.txt** file.

## Retrieving gene orthologs from mammalian genomes
For getting orthologues from a list of mammalian genomes, I used the BiomaRt library implemented in R. The R script **get_orthologues_biomart.R** was written to perform the retrieving task.

