###Build bs_data
###Alberto Vicens

library(magrittr)
setwd("/home/uvi/be/avs/lustre/cancer_genes_selection")

#1. Load file with species codes
spp<-read.table("data/species_codes.txt",sep=" ",header=FALSE)
spp<-spp[,-2]
names(spp) <- c("fullName","SppCode")
spp$SppCode<-as.character(spp$SppCode)
spp$SppCode<-gsub("ENSP00","ENSG00",spp$SppCode)

#Create dataframe which LRTs will be added to
bs<-data.frame(Species=spp$fullName)

#2. Load file with node IDs
ids<-read.table("genes/CCND2/paml/CCND2_node_ids_parsed.txt", sep=";",col.names = "GeneID")
ids$nodeID<-rownames(ids)
ids$SppCode<-substr(ids$GeneID,1,6)

#3. Load LRT file
lrt<-read.table("genes/CCND2/paml/CCND2_LRT.txt", sep="|", header=TRUE)
lrt$nodeID<- lrt$Alternativemodel 
lrt$nodeID %<>% gsub("1$","",.) %>% gsub("[^0-9]+","",.)

#4. Merge dataframes by SppCode
merged1<-merge(ids[,c("SppCode","nodeID")],lrt[,c("nodeID","p.value")],by="nodeID")
merged2<-merge(spp,merged1,by="SppCode", all.x=TRUE)
merged2$p.value<-as.numeric(gsub("\\*","",merged2$p.value))
merged2$p.adj<-p.adjust(p=merged2$p.value,method = "fdr")

#5. Add LRT data to the branch-site dataframe
bs<-cbind(bs,merged2$p.adj)
colnames(bs)[ncol(bs)]=genename

