# Molecular Adaption of Cancer Genes
## Description
In this repository I describe the workflow followed to process and analyze the data of the study "Selective constraints driving the Germinal Evolution of Cancer Genes in Mammals" *(manuscript in preparation)*.
The repository contains the scripts and some example data used in the study.

## Cancer genes collection
I retrieved cancer genes from the Cancer Gene Census (CGC) repository of the COSMIC database. A list of 574 genes classified as Tier 1 (i.e. those genes with a documented activity in cancer) was downloaded. This dataset can be accessed in the link https://cancer.sanger.ac.uk/cosmic/census?tier=1 and was downloaded on March 5th, 2018.
The dataset was cleaned in order to be used for later analysis with the R script *build_fundata_table.R* 

The variables included in the dataset are:
* Somatic
* Germline
* Tumor type
* Tissue type
* Molecular genetics
* Role in cancer
* Mutation type