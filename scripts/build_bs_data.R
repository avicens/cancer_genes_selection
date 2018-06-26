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

#2. Get significance values of LRTs for each gene

##2.1 Get a list with genenames
genedir<-"/home/uvi/be/avs/lustre/cancer_genes_selection/genes"
genelist<-list.files(genedir)

##2.2 Run a loop to extract p-adjusted values of each gene and attach to bs dataframe

for (i in 2) {
###2.2.1 Load file with node IDs
  ids<-read.table(file=paste(genedir,"/",genelist[i],"/paml/",genelist[i],"_nodeIDs_parsed.txt", 
                      sep=""),sep=":",col.names = "GeneID")
  ids$nodeID<-rownames(ids)
  ids$SppCode<-substr(ids$GeneID,1,6)


###2.2.2. Load LRT file
lrt<-read.table(file=paste(genedir,"/",genelist[i],"/paml/",genelist[i],"_LRT.txt", sep=""), sep="|", header=TRUE)
lrt$nodeID<- lrt$Alternativemodel 
lrt$nodeID %<>% gsub("1$","",.) %>% gsub("[^0-9]+","",.)

#4. Merge dataframes by SppCode
merged1<-merge(ids[,c("SppCode","nodeID")],lrt[,c("nodeID","p.value")],by="nodeID")
merged2<-merge(spp,merged1,by="SppCode", all.x=TRUE)
merged2$p.value<-as.numeric(gsub("\\*","",merged2$p.value))
merged2$p.adj<-p.adjust(p=merged2$p.value,method = "fdr")

#5. Add LRT data to the branch-site dataframe
bs<-cbind(bs,merged2$p.adj)
colnames(bs)[ncol(bs)]<-genelist[i]
}

