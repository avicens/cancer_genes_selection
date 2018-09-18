#Open data frame with COSMIC information
setwd("~/Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/")

#1. Load table downloaded from COSMIC
fundata<-read.csv("~/Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/tables/cosmic_cancer_genes.csv", header=TRUE)

#2. Remove variables which will not be analyzed (I)
fundata<-fundata[,-c(2,3,5:7,10:12,17:20)]
colnames(fundata)<-c("gene","location", "somatic","germline","tissue.type","inheritance","cancer.role","impact")

#3. Curate mutation columns
fundata$somatic<-sub("^$","no",fundata$somatic)
fundata$germline<-sub("^$","no",c(as.character(fundata$germline)))

mut<-character()
for (i in 1:nrow(fundata)) {
  
  if (fundata$germline[i]=="yes" & fundata$somatic[i]=="no") {mut[i]<-"germinal"}
  if (fundata$germline[i]=="no" & fundata$somatic[i]=="yes") {mut[i]<-"somatic"}  
  if (fundata$germline[i]=="yes" & fundata$somatic[i]=="yes") {mut[i]<-"som&germ"}
}
fundata$mut.class<-as.factor(mut)

#4. Extract chromosome number and chromosome type
chr<-sapply(strsplit(as.character(fundata$location),split=":"),"[",1); chr<-gsub("\"X,Y","X",chr,ignore.case = T)
fundata$chr<-chr

chrtype=chr; id<-which(chrtype!="X"); chrtype[id]="Autosome"
fundata$chr.type=chrtype

#5. Curate inheritance column
fundata$inheritance<-gsub("Dom", "Dominant", fundata$inheritance) 
fundata$inheritance<- gsub("Rec", "Recessive", fundata$inheritance)
fundata$inheritance<-as.factor(fundata$inheritance)

#6. Remove variables which will not be analyzed (II)
fundata<-fundata[,-c(3,4,10)]

