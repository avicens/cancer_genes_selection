#Open data frame with COSMIC information
fundata<-read.table("Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/tables/cosmic_cancer_genes_compacted.tsv",sep="\t",header=T)
colnames(fundata)<-c("gene","name","somatic","germline","tumour.type","tissue.type","genetics","cancer.role","mutation.type")
                                                   
#Curate mutation columns
fundata$somatic<-sub("^$","no",fundata$somatic)
fundata$germline<-sub("^$","no",c(as.character(fundata$germline)))

mut<-character()
for (i in 1:nrow(fundata)) {
  
  if (fundata$germline[i]=="yes" & fundata$somatic[i]=="no") {mut[i]<-"germinal"}
  if (fundata$germline[i]=="no" & fundata$somatic[i]=="yes") {mut[i]<-"somatic"}  
  if (fundata$germline[i]=="yes" & fundata$somatic[i]=="yes") {mut[i]<-"som&germ"}
  }
fundata$mut.class<-as.factor(mut)

#Curate tissue.type column
fundata$tissue.type<-as.factor(gsub("\"| ","",fundata$tissue.type))
tissues<-strsplit(levels(fundata$tissue.type),split=",|;")
tissues<-unique(unlist(tissues))

tt<-unlist(sapply(as.character(tissue.type),strcomb),use.names=F)
fundata$tissue.type<-as.factor(sub(" ","",tt))

#Curate tumor.type column
fundata$tumour.type<-as.factor(gsub("\"","",fundata$tumour.type))

#Add data from MartinCorena
mcdata<-read.table("Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/data/datos_iñigo.tsv",sep="\t",header=T)
mcdata<-mcdata[c(1:202),-ncol(mcdata)]
driver_mc<-sapply(gene,function(x) {x %in% mcdata$gene})

mcdnds<-read.table("Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/data/dNdScv_output_PANCANCER.txt",sep="\t",header=T)
idy <- sapply(as.character(fundata[,1]),function (x) which(as.character(mcdnds[,1])==x))
mcdnds2<-mcdnds[idy,]
fundata$dnds_mc<-mcdnds2$wmis3