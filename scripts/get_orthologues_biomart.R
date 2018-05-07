library(biomaRt)

#Load database
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Obtain list of orthologues
spfile<-read.table("/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection/data/species_ensembl_code.txt",header = F)
splist<-paste(spfile$V1,"_homolog_ensembl_peptide",sep="")
splist=splist[-11] #Discard human

#Obtain list of human peptide ID
pepfile<-read.table("/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection/data/protein_ensembl_id.txt",header = F)
peplist=pepfile$V1

#Set attributes
attributes = c("external_gene_name","ensembl_peptide_id")
orthologues = getBM(attributes,filters="ensembl_peptide_id",values=peplist,TRUE, mart = human, uniqueRows=TRUE)
colnames(orthologues)=c("Gene","hsapiens")

for (i in 1:length(splist)) {
ortho=getBM(c("external_gene_name",splist[i]),filters="ensembl_peptide_id",values=peplist, mart = human, uniqueRows=TRUE)
ortho=ortho[!duplicated(ortho$external_gene_name),2]
orthologues=cbind(orthologues,ortho)
spname<-sapply(strsplit(splist[i],"_"), "[",1)
colnames(orthologues)[ncol(orthologues)]=spname
}

#Extract ortholog lists
for (j in 1:nrow(orthologues)){
gene=orthologues[j,1]
ortholist=orthologues[j,-1]
ortholist=ortholist[!ortholist[1,]==""]
write.table(ortholist,file=paste("/mnt/lustre/scratch/home/uvi/be/avs/cancer_genes_selection/genes//",gene,"/",gene,"_orthologs2.txt",sep=""),quote = F,col.names = F, row.names = F, sep="\n")
}
