library(biomaRt)

##Polymorphism data
#Load database
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pepfile<-read.table("/home/user/Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/data/protein_ensembl_id.txt",header = F)
peplist=pepfile$V1

#Set attributes
filt = c("germ_line_variation_source","ensembl_peptide_id", "so_mini_parent_name")
val=list(glvs="dbSNP", epi=peplist, smpn = "coding_sequence_variant")
attr = c("external_gene_name","variation_name","synonymous_status")
polim<-getBM(attributes=attr,filters= filt, values = val, mart =human, uniqueRows = T)

#Split the table for genes
s<-split(polim, polim$ensembl_gene_id)

#Get number of synonymous (S) and nonsynonymous (N) polymorphims
N<-sapply(s,function(x) length(grep("missense",x[,3])))
S<-sapply(s,function(x) length(grep("synonymous",x[,3])))
polimNS<-cbind(N,S)

##Divergence data
filt2 = "ensembl_peptide_id"
attr2=c("external_gene_name","ptroglodytes_homolog_dn","ptroglodytes_homolog_ds")
div<-getBM(attributes=attr2,filters= filt2, values = peplist, mart=human, uniqueRows = T)
