# Positive Selection of Cancer Genes
## Description
In this repository I describe the workflow followed to process and analyze the data of the study "Selective pressures on human cancer genes along the evolution of mammals" (Genes 2018, 9(12), 582; https://doi.org/10.3390/genes9120582).

## Repository content
The repository contains the scripts written to perform the computational tasks of this study, as well as the input files required as input in some steps.

The scripts were written to upload jobs to a SLURM-based scheduling system (CESGA Finis Terrae II). Those scripts whose filename end with *array.sh* call to launch parallel tasks.

## Sequence Data Collection
### Cancer genes dataset
I retrieved cancer genes from the [Cancer Gene Census Database](https://cancer.sanger.ac.uk/cosmic/census?tier=1) of the COSMIC v84 database. A list of **574 genes** classified as Tier 1 (those genes with a documented activity in cancer) was downloaded. 
This dataset was accessed on March 5th, 2018, and corresponds to the **Table S1** of the manuscript. 

The dataset was cleaned in order to be used for later analysis with the R script **build_fundata_table.R** 

The variables included in the dataset are:

* Mutation type (Somatic/Germline)
* Inheritance (Dominant/Recessive)
* Cancer role 
* Tissue type


### Retrieving human protein annotations
The first column of the cancer genes dataset, which contains the gene names, was used to retrieve information for human genes from [Ensembl BioMArt](http://www.ensembl.org/biomart/martview/4ee102879139fc3bf745f9a867064956).

In the exported file, genes with no protein annotated were discarded, obtaining **535 genes**. On these genes, the protein isoform with the best transcript support level was chosen.

```{bash, filter_genes}
tail -n +2 ensembl_gene_list.tsv| grep "\<tsl1\>" | grep 'GENCODE'| sort -k4,4 -k6nr,6 | sort -u -k 4,4 > ensembl_gene_uniq.tsv
```

The obtained dataset is the file **ensembl_gene_uniq.tsv**.

In this step, and in order to apply posterior analysis for each gene, a working folder was created for each gene with the bash script **create_gene_folders.sh**.

The Ensembl_protein_ID column was extracted from the **ensembl_gene_uniq.tsv** table, getting the **protein_ensembl_id.txt** file.
```{bash}
cut -f2 ensembl_gene_uniq.tsv > protein_ensembl_id.txt
```
### Collecting gene orthologs from mammalian genomes
For getting orthologues from a list of 32 mammalian genomes (see **Table S3** of the manuscript), I used the BiomaRt library implemented in R. The R script **get_orthologues_biomart.R** was written to perform this task.

Ortholog coding sequences from all species were downloaded with the [BioMart Perl API](http://www.ensembl.org/info/data/biomart/biomart_perl_api.html#biomartperl).

For this task, I wrote the **download_ortholog_cds_array.sh** script that call to the **download_cds_from_protid.pl** perl script iteratively for each gene and species.

## Multiple sequence alignment (MSA)
Coding nucleotide sequences were aligned using the software MACSE (*Ranwez et al. 2011*). This program accounts for frameshifts and stop codons, and it is optimal for aligning coding sequences. To perform the MSA for every gene, I wrote the script **multiple_alignment_array.sh**

### Curating MSA
In the MSA curation process, the following steps were applied:

* Remove STOP codons
* Clean up secondary annotations of fasta headers to jusy retain Ensembl Gene IDs
* Discarding duplicate sequences
* Substitution of "!" characters (inserted by MACSE to corect frameshits) for "N" characters.

Once the MSAa were curated, they were trimmed to remove poorly aligned positions and sequences using TrimAl (*Capella-Gutierrez et al. 2009*).
The applied  parameters were:

* Delete columns with gaps in more than 60% of the sequences.
* Delete columns with a similarity score lower than 0.1
* Remove sequences not covering at least the 60% of residues that achieve an overlap, with the rest of the sequences, of 0.75.

All these steps were sequentially performed with the script **process_msa.sh**

## Tree reconstruction
Phylogenetic trees for each gene were built using the program RAxML-ng (*Stamatakis 2014, Kozlov 2016*). The phylogenetic reconstruction analysis included ML tree search + non-parametric bootstrap, with these parameters:

  * 10 randomized parsimony starting trees
  * General Time Reversible substitution model with discrete GAMMA model of rate heterogeneity with 4 categories (GTR+G)
  * 100 bootstrap replicates
  
The code to build phylogenetic trees for every gene was written into **run_raxml_array.sh**

### Estimation of tree distances
In order to retain reliable gene families, I estimated the distance of each gene tree with respect to a well supported mammalian species tree using ETE3 (*Huerta-Cepas et al. 2016*). The code to launch this analysis is written in **compare_trees.sh**

## Selection analysis
Tests of positive selection were performed with the Codeml program implemented in PAML (Yang 2007), but I used the ETE3 evol package (Huerta-Cepas et al. 2006) to conduct analysis. The set parameters in common for all models were:

  * Codon frequency: codon table
  * Initial omega: 0.7
  * Number of gamma categories: 4

Using MSA and gene trees as input files, I applied the following selection site models on each gene:

  * Estimation of global evolutionary rate (dN/dS) with the model **M0**.
  * Test of variation in evolutionary rate across sites with the model **M3**
  * Test of positive slection applying models **M1 vs M2** and **M8 vs M8a**

All these models were conducted with the script **run_codeml.sh**
