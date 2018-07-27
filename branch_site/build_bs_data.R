###Build bs_data
###Alberto Vicens

library(magrittr)
setwd("/home/uvi/be/avs/lustre/cancer_genes_selection/")

#1. Load file with species codes
spp<-read.table("data/species_codes.txt",sep=" ",header=FALSE)
spp<-spp[,-2]
names(spp) <- c("fullName","SppCode")
spp$SppCode<-as.character(spp$SppCode)
spp$SppCode<-gsub("ENSP00","ENSG00",spp$SppCode)

#Create dataframe which LRTs will be added to
bs<-data.frame(Species=spp$fullName)

#2. Get significance values of LRTs for each gene

##2.1 Get a list with gene names
genedir<-"/home/uvi/be/avs/lustre/cancer_genes_selection/genes"
genelist<-list.files(genedir)

##2.2 Run a loop to extract p-adjusted values of each gene and attach them to the bs dataframe

for (i in 1:length(genelist)) {

  ###2.2.1 Load file with node IDs
  ids<-read.table(file=paste(genedir,"/",genelist[i],"/paml/",genelist[i],"_nodeIDs_parsed.txt", 
                      sep=""),sep=":",col.names = "GeneID")
  ids$nodeID<-rownames(ids)
  ids$SppCode<-substr(ids$GeneID,1,6)


  ###2.2.2. Load LRT file
  lrt<-read.table(file=paste(genedir,"/",genelist[i],"/paml/",genelist[i],"_LRT.txt", sep=""), sep="|", header=TRUE)
  lrt$nodeID<- lrt$Alternativemodel 
  lrt$nodeID %<>% gsub("1$","",.) %>% gsub("[^0-9]+","",.)

  #2.2.3. Merge dataframes by SppCode
  merged1<-merge(ids[,c("SppCode","nodeID")],lrt[,c("nodeID","p.value")],by="nodeID")
  merged2<-merge(spp,merged1,by="SppCode", all.x=TRUE)
  
  #2.2.4. Estimate p-value and p-adj (multiple testing)
  merged2$p.value<-as.numeric(gsub("\\*","",merged2$p.value))
  merged2$p.adj<-p.adjust(p=merged2$p.value,method = "fdr")

  #2.2.5. Add LRT data to the branch-site dataframe
  bs<-cbind(bs,merged2$p.adj)
  colnames(bs)[ncol(bs)]<-genelist[i]
}

#Export table with p-adj per gene/species
write.table(bs, file="branch_site/bs_p-adj.tsv", sep="\t", quote = F, row.names = F, col.names = T)

#Number of analyzed genes per species
ngenes<-apply(bs,1, function(x) length(x[!is.na(x[-1])]))

#Extract list of positively selected genes per species
psglist<-list()

for (j in 1:nrow(bs)) {
r<-bs[j,]; r<-r[,!is.na(r)]
psglist<-c(list(ps=names(r[,which(r< 0.01)])),psglist)
}

names(psglist)<-bs[,1]

#Number of positively selected genes per species
psgenes<-unlist(lapply(psglist,function(x) length(x)))
percgenes<-(psgenes/ngenes)*100

ps<-as.data.frame(rbind(psgenes,ngenes))
ps<-ps[idx]
ps$genes<-c("PS", "Total")
ps.long<-melt(ps,id.vars = "genes")
names(ps.long)<-c("GeneClass","Species","Count")

gpl<-ggplot(ps.long,aes(x=Species, y=Count, fill=GeneClass)) 
gpl+ geom_bar(stat="identity", position = "dodge") + coord_flip() +
  scale_fill_discrete(name="Genes", labels=c("Positive selected","Total"))



