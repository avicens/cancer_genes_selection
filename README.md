# Molecular Adaption of Cancer Genes
## Description
In this repository I describe the workflow followed to process and analyze the data of the study "Selective constraints driving the Germinal Evolution of Cancer Genes in Mammals" *(manuscript in preparation)*.

## Repository content
The repository contains the scripts written to perform the computational tasks of this study, as well as the input files required as input in some steps.

The scripts were written to upload jobs to a SLURM-based scheduling system (CESGA Finis Terrae II). Those scripts whose filename end with *array.sh* call to launch parallel tasks.

## Sequence Data Collection
### Cancer genes dataset
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

### Retrieving human protein annotations
The first column of the cancer genes dataset, which contains the gene names, was extracted and the information for human genes was retrieved using Ensembl BioMArt (http://www.ensembl.org/biomart/martview/4ee102879139fc3bf745f9a867064956). In the exported file, genes with no protein annotated were discarded, obtaining 535 genes. On these genes, the protein isoform with the best transcript support level was chosen. The obtained dataset is the file **ensembl_gene_uniq.tsv**.

In this step, and in order to apply the later analysis per gene, a working folder was created for each gene with the bash script **create_gene_folders.sh**.

The Ensembl_protein_ID column was extracted from the **ensembl_gene_uniq.tsv** table, getting the **protein_ensembl_id.txt** file.

### Collecting gene orthologs from mammalian genomes
For getting orthologues from a list of 32 mammalian genomes (see Table S2 of the manuscript), I used the BiomaRt library implemented in R. The R script **get_orthologues_biomart.R** was written to perform this task.
Ortholog coding sequences from all species for each gene were downloaded with BioMart Perl API. For this task, I wrote the **download_ortholog_cds_array.sh** script that call to the **download_cds_from_protid.pl** perl script iteratively for each gene and species.

## Multiple sequence alignment (MSA)
Coding nucleotide sequences were aligned using the software MACSE (*Ranwez et al. 2011*). This program accounts for frameshifts and stop codons, and it is optimal for aligning coding sequences. To perform the MSA for every gene, I wrote the script **multiple_alignment_array.sh**

### Curating MSA
In the MSA curation process, were applied the following steps:

* Clean up secondary annotations of fasta headers to jusy retain Ensembl Gene IDs
* Discarding duplicate sequences
* Substitution of "!" characters (inserted by MACSE to corect frameshits) for "N" characters.

Every step were sequentially performed with the script **curate_msa.sh**

### Trimming MSA
Once the MSAa were curated, they were trimmed to remove poorly aligned positions and sequences using TrimAl (*Capella-Gutierrez et al. 2009*).
The applied  parameters were:

* Delete columns with gaps in more than 60% of the sequences.
* Delete columns with a similarity score lower than 0.1
* Remove sequences not covering at least the 60% of residues that achieve an overlap, with the rest of the sequences, of 0.75.

The trimming task was implemented in the script **trim_msa.sh** 

## Tree reconstruction
Phylogenetic trees for each gene were built using the program RAxML-ng (*Stamatakis 2014, Kozlov 2016*). The phylogenetic reconstruction analysis included ML tree search + non-parametric bootstrap, with these parameters:

  * 10 randomized parsimony starting trees
  * General Time Reversible substitution model with discrete GAMMA model of rate heterogeneity with 4 categories (GTR+G)
  * 100 bootstrap replicates

