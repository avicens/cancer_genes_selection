#Combine functional and evolutionary data frames
setwd("Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/") #Set working directory
load("data/evoldata") #Load data frame with evolutionary data (PAML)
load("data/fundata") #Load data frame with functional data (COSMIC)

fulldata<-merge(evoldata, fundata,by="gene")
attach(fulldata)

#load functionss and libraries
source("scripts/dndspercat.R") #Load function to generate dataframes per category
source("scripts/make_comparisons.R") #Load function to generate comparisons invariables with more than 2 classes
library(ggpubr) #Load "ggpubr" library for getting labelled box plot

##Casting tables

#Mutations type (somatic, germline)
mutclasses<-unique(unlist(strsplit(as.character(fulldata$mut.class),split=", ")))
mutclassnames<-c("Somatic","Som + Germ", "Germinal")
mutclassdf<-dndspercat(mutclasses,mutclassnames,fulldata,mut.class)
mutclasscomp<-make.comparisons(mutclassnames) #Generate a list with all possible comparisons

#Genetics
genet<-levels(as.factor(fulldata$inheritance))
genet<-c("Dom","Rec")
genames<-c("Dominant","Recessive")
genout<-dndspercat(genet,genames,fulldata,inheritance)

#Cancer role
canrole<-gsub("\"","",levels(fulldata$cancer.role))
canrole<-unique(unlist(strsplit(canrole,", ")))
canout<-dndspercat(canrole, canrole,fulldata,cancer.role)

#Mutation type (functional)
mutype<-gsub("\"","",fulldata$impact)
mutype<-unique(unlist(strsplit(mutype,", |. |,| ")))
mutype2<-c("Mis","N")
mutnames<-c("Missense","Nosense")
mutout<-dndspercat(mutype,mutnames,fulldata,impact)

##Tissue type
tissues<-unique(unlist(strsplit(levels(fulldata$tissue.type),split=",|;")))
tisnames<-c("Epithelial","Leukemia/Lymphoma","Mesenchimal","Other")
tisdf<-dndspercat(tissues,tisnames,fulldata,tissue.type)

##5.Cancer role
canrole<-gsub("\"","",levels(fulldata$cancer.role))
canrole<-unique(unlist(strsplit(canrole,", ")))
candf<-dndspercat(canrole, canrole,fulldata,cancer.role)

###Statistical comparisons
ger_plot<-ggboxplot(fulldata, x = "germline", y = "avomega",color="germline",palette="jco",ylab = "Global dN/dS") + stat_compare_means(method="wilcox.test",label.x=1.5,cex=5)
ger_plot2<-ggboxplot(fulldata, x = "germline", y = "avomega",color="germline",palette="jco",facet.by="somatic",ylab = "Global dN/dS", short.panel.labs = F) + stat_compare_means(method="wilcox.test",label.x=1.5,cex=4)

ggarrange(ger_plot,ger_plot2, ncol=2,nrow=1,widths=c(2,3))


mosaicplot(table(fulldata$mut.class,fulldata$PS),ylab="Positive selection",color = c("lightblue","red"),main = "Mutation type",cex.axis =1,xlab=paste("p-value(chisq.test)",round(chisq.test(table(fulldata$mut.class,fulldata$PS))$p.value,3),sep="="))

mosaicplot(table(genout$class,genout$PS),ylab="Positive selection",color = c("lightblue","red"),main = "Inheritance",cex.axis =1,xlab=paste("p-value(chisq.test)",round(chisq.test(table(genout$class,genout$PS))$p.value,3),sep="="))

mosaicplot(table(mutout$class,mutout$PS),ylab="Positive selection",color = c("lightblue","red"),main = "Functional impact",cex.axis =1,xlab=paste("p-value(chisq.test)",round(chisq.test(table(mutout$class,mutout$PS))$p.value,3),sep="="))

mosaicplot(table(tisout$class,tisout$PS),ylab="Positive selection",color = c("lightblue","red"),main = "Tissue type",cex.axis =1,xlab=paste("p-value(chisq.test)",round(chisq.test(table(tisout$class,tisout$PS))$p.value,3),sep="="))

mosaicplot(table(canout$class,canout$PS),ylab="Positive selection",color = c("lightblue","red"),main = "Cancer role",cex.axis =1,xlab=paste("p-value(chisq.test)",round(chisq.test(table(canout$class,canout$PS))$p.value,3),sep="="))

mosaicplot(table(fulldata$chr.type,fulldata$PS),ylab="Positive selection",color = c("lightblue","red"),main = "Chromosome type",cex.axis =1,xlab=paste("p-value(chisq.test)",round(chisq.test(table(fulldata$chr.type,fulldata$PS))$p.value,3),sep="="))
```

We also compared evolution of cancer genes between somatic and germinal levels. For this we took data from Martín-Corena et al. (2017) and tested wheter global dN/dS vary in genes under positive selection at the somatic level.
```{r ps_plot, eval=T, fig.align="center", dpi=100, fig.width=7, fig.height=14,fig.cap="Figure . Up: Boxplot comparing global dN/dS for genes previously identified (or not) under positive somatic selection. Statistical significance of comparison is indicated. Botton: correlation between germinal and somatic (missense) dN/dS. Blue points represent genes under soamtic positive selection according to MartínCorena et al. 2017." }
attach(fulldata)
psplot<-ggboxplot(fulldata, x = "driver_mc", y = "avomega",color="driver_mc",palette="jco",ylab = "Global dN/dS", short.panel.labs = F,xlab="Positive selection",submain="Driver genes from Martín-Corena et al. 2017") + stat_compare_means(method="wilcox.test",label.x=1.5,cex=5)

somgerplot<-ggplot(fulldata,aes(avomega,log10(dnds_mc),colour=driver_mc))+geom_point() + xlab("germinal dN/dS") +ylab("log10(somatic dN/dS)") + annotate(x=0.5, y=1.5, label=paste("R = ", round(cor(avomega, log10(dnds_mc)),2)),geom="text", size=5)

ggarrange(psplot,somgerplot,nrow=2,ncol=1)
```

