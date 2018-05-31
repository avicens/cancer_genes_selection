#Combine functional and evolutionary data frames
setwd("Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/") #Set working directory
load("data/evoldata") #Load data frame with evolutionary data (PAML)
load("data/fundata") #Load data frame with functional data (COSMIC)

fulldata<-merge(evoldata, fundata,by="gene")
attach(fulldata)


#Boxplot
#dN/dS
driv_dNdS<-ggboxplot(fulldata, x = "driver", y = "dNdS",color="driver",palette="jco",ylab = "Global dN/dS", xlab=FALSE, legend=NULL) 
driv_dNdS2<- driv_dNdS + stat_compare_means(method="t.test",label.x=1.5,cex=5)

#dN
driv_dN<-ggboxplot(fulldata, x = "driver_mc", y = "dN",color="driver_mc",palette="jco",ylab = "dN", xlab=FALSE, legend=NULL) 
driv_dN2<- driv_dN + stat_compare_means(method="t.test",label.x=1.5,cex=5)

#dS
driv_dS<-ggboxplot(fulldata, x = "driver_mc", y = "dS",color="driver_mc",palette="jco",ylab = "dS", xlab=FALSE, legend=NULL) 
driv_dS2<-driv_dS + stat_compare_means(method="t.test",label.x=1.5,cex=5)


### Comparison of dN/dS across categories
source("scripts/dndspercat.R") #Load function to generate dataframes per category
source("scripts/make_comparisons.R") #Load function to generate comparisons invariables with more than 2 classes

library(ggpubr) #Load "ggpubr" library for getting labelled box plot

###1.Molecular constraints

##1.1. Mutation type (somatic, germline) I: comparison of a variable after grouping its data with the oter variable
#Comparing somatic mutations
som_plot<-ggboxplot(fulldata, x = "somatic", y = "avomega",color="somatic",palette="jco",ylab = "Global dN/dS", short.panel.labs = F) + stat_compare_means(method="wilcox.test",label.x=1.5,cex=5)
somfacetpl<-ggboxplot(fulldata, x = "somatic", y = "avomega",color="somatic",palette="jco",facet.by="germline",ylab = "Global dN/dS", short.panel.labs = F) + stat_compare_means(method="wilcox.test",label.x=1.5,cex=4)
#Comparing germinal mutations
ger_plot<-ggboxplot(fulldata, x = "germline", y = "avomega",color="germline",palette="jco",ylab = "Global dN/dS") + stat_compare_means(method="wilcox.test",label.x=1.5,cex=5)
gerfacetpl<-ggboxplot(fulldata, x = "germline", y = "avomega",color="germline",palette="jco",facet.by="somatic",ylab = "Global dN/dS", short.panel.labs = F) + stat_compare_means(method="wilcox.test",label.x=1.5,cex=4)

#Mutations type (somatic, germline) II: comparison of 3 groups
mutclasses<-unique(unlist(strsplit(as.character(fulldata$mut.class),split=", ")))
mutclassnames<-c("Somatic","Som + Germ", "Germinal")
mutclassdf<-dndspercat(mutclasses,mutclassnames,fulldata,mut.class)

#Data exploration
table(mutclassdf$class) #Check the sample size for each group
mutclasscomp<-make.comparisons(mutclassnames) #Generate a list with all possible comparisons

#Boxplot
#dN/dS
mutclass_dNdS<-ggboxplot(mutclassdf, x = "class", y = "dNdS", select= mutclassnames, color="class",palette="jco",xlab=FALSE, ylab = "Global dN/dS") 
mutclass_dNdS2<-mutclass_dNdS + stat_compare_means(method = "t.test", comparisons = mutclasscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 1, label.x = 0.7)

#dN
mutclass_dN<-ggboxplot(mutclassdf, x = "class", y = "dN", select= mutclassnames, color="class",palette="jco",xlab=FALSE, ylab = "dN") 
mutclass_dN2<-mutclass_dN + stat_compare_means(method = "t.test", comparisons = mutclasscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 4, label.x = 0.7)

#dS
mutclass_dS<-ggboxplot(mutclassdf, x = "class", y = "dS", select= mutclassnames, color="class",palette="jco",xlab=FALSE, ylab = "dS") 
mutclass_dS2<-mutclass_dS + stat_compare_means(method = "t.test", comparisons = mutclasscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 13, label.x = 0.7)


##1.2. Molecular genetics
genet<-levels(fulldata$genetics)
genet<-c("Dom","Rec")
genames<-c("Dominant","Recessive")
gendf<-dndspercat(genet,genames,fulldata,genetics)

#Data exploration
table(gendf$class) #Check the sample size for each group
gencomp<-make.comparisons(genames) #Generate a list with all possible comparisons

#dN/dS
gen_dnds<-ggboxplot(gendf, x = "class", y = "dNdS", select= genames, color="class",palette="jco",xlab=FALSE, ylab = "Global dN/dS") 
gen_dnds2<-gen_dnds + stat_compare_means(method = "t.test", label.y = 1, label.x = 0.7, cex=5)

#dN
gen_dN<-ggboxplot(gendf, x = "class", y = "dN", select= genames, color="class",palette="jco",xlab=FALSE, ylab = "dN") 
gen_dN2<-gen_dN + stat_compare_means(method = "t.test", label.y = 5, label.x = 0.7)

#dS
gen_dS<-ggboxplot(gendf, x = "class", y = "dS", select= genames, color="class",palette="jco",xlab=FALSE, ylab = "dS") 
gen_dS2<-gen_dS + stat_compare_means(method = "t.test", label.y = 15, label.x = 0.7)

#1.3. Functional effect of coding mutations
mutype<-gsub("\"","",fulldata$mutation.type)
mutype<-unique(unlist(strsplit(mutype,", |. |,| ")))
mutype2<-c("Mis","N")
mutnames<-c("Missense","Nonsense")
mutefdf<-dndspercat(mutype2,mutnames,fulldata,mutation.type)

#Data exploration
table(mutefdf$class) #Check the sample size for each group

#Boxplot
#dN/dS
mutef_dnds<-ggboxplot(mutefdf, x = "class", y = "dNdS", select= mutnames, color="class",palette="jco",xlab=FALSE, ylab = "Global dN/dS") 
mutef_dnds2<- mutef_dnds + stat_compare_means(method = "t.test", label.y = 0.6, label.x = 1.5, cex= 5)

#dN
mutef_dN<-ggboxplot(mutefdf, x = "class", y = "dN", select= mutnames, color="class",palette="jco",xlab=FALSE, ylab = "dN") 
mutef_dN2<-mutef_dN + stat_compare_means(method = "t.test", label.x= 0.5)

#dS
mutef_dS<-ggboxplot(mutefdf, x = "class", y = "dS", select= mutnames, color="class",palette="jco",xlab=FALSE, ylab = "dS") 
mutef_dS2<- mutef_dS + stat_compare_means(method = "t.test", label.x = 0.5)


###2. Cancer constraints###
###2.1. Tissue type
tissues<-unique(unlist(strsplit(levels(fulldata$tissue.type),split=",|;")))
tisnames<-c("Epithelial","Leuk/Lymp","Mesen","Other")
tisdf<-dndspercat(tissues,tisnames,fulldata,tissue.type)

#Data exploration
table(tisdf$class) #Check the sample size for each group
tiscomp<-make.comparisons(tisnames) #Generate a list with all possible comparisons

#Boxplot
#dN/dS
tis_dnds<-ggboxplot(tisdf, x = "class", y = "dNdS", select= tisnames, color="class",palette="jco",xlab=FALSE, ylab = "Global dN/dS") 
tis_dnds2<-tis_dnds + stat_compare_means(method = "t.test", comparisons = tiscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 1.2, label.x = 0.7)

#dN
tis_dN<-ggboxplot(tisdf, x = "class", y = "dN", select= tisnames, color="class",palette="jco",xlab=FALSE, ylab = "dN") 
tis_dN2<-tis_dN + stat_compare_means(method = "t.test", comparisons = tiscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 5, label.x = 0.7)

#dS
tis_dS<-ggboxplot(tisdf, x = "class", y = "dS", select= tisnames, color="class",palette="jco",xlab=FALSE, ylab = "dS") 
tis_dS2<-tis_dS + stat_compare_means(method = "t.test", comparisons = tiscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 15, label.x = 0.7)

##2.2.Cancer role
canrole<-gsub("\"","",levels(fulldata$cancer.role))
canrole<-unique(unlist(strsplit(canrole,", ")))
candf<-dndspercat(canrole, canrole,fulldata,cancer.role)

#Data exploration
table(candf$class) #Check the sample size for each group
cancomp<-make.comparisons(canrole) #Generate a list with all possible comparisons

#dN/dS
can_dnds<-ggboxplot(candf, x = "class", y = "dNdS", select= canrole, color="class",palette="jco",xlab=FALSE, ylab = "Global dN/dS") 
can_dnds2<-can_dnds + stat_compare_means(method = "t.test", comparisons = cancomp, label ="p.signif") + stat_compare_means(method = "anova", label.y = 1, label.x = 0.7, cex=5)

#dN
can_dN<-ggboxplot(candf, x = "class", y = "dN", select= canrole, color="class",palette="jco",xlab=FALSE, ylab = "dN") 
can_dN2<- can_dN + stat_compare_means(method = "t.test", comparisons=cancomp, label ="p.signif") + stat_compare_means(method = "anova", label.y = 4.2, label.x = 0.7) 

#dS
can_dS<-ggboxplot(candf, x = "class", y = "dS", select= canrole, color="class",palette="jco",xlab=FALSE, ylab = "dS")
can_dS2<-can_dS + stat_compare_means(method = "t.test", comparisons = cancomp, label ="p.signif") + stat_compare_means(method = "anova", label.y = 4.2, label.x = 0.7)

###2.3. Tumour type
tumour<-unique(unlist(strsplit(levels(fulldata$tumour.type),split=",|;")))
tumger<-unique(unlist(strsplit(levels(fundata2$Tumour.Types.Germline.),split=",|;")))
tumger<-tumger[!(tumger=="FANCD1)")]


##multiplot (dN/dS)
ggarrange(mutclass_dNdS2, gen_dnds2, mutef_dnds2, tis_dnds2, can_dnds2, driv_dNdS2, nrow=2, ncol=3)

##Multiplot (dN,dS)
ggarrange(mutclass_dN2,mutclass_dS2,gen_dN2,gen_dS2,mutef_dN2, mutef_dS2,tis_dN2,tis_dS2,can_dN2,can_dS2, driv_dN2, driv_dS2, nrow=2, ncol= 6)

#Comparison of proportion of positively selected genes
mutclassChisq<-chisq.test(table(mutclassdf$psgenes, mutclassdf$class))
genChisq<-chisq.test(table(gendf$psgenes, gendf$class))
mutefChisq<-chisq.test(table(mutefdf$psgenes, mutefdf$class))
tisChisq<-chisq.test(table(tisdf$psgenes, tisdf$class))
canChisq<-chisq.test(table(candf$psgenes, candf$class))
drivChisq<-chisq.test(table(fulldata$psgenes,fulldata$driver))

#Mosaic plost comparing proportion of positively selected genes
png(file="figures/mosaicplots1")
par(mfrow=c(2,3), oma=c(rep(1,4)))
mosaicplot(table(mutclassdf$class,mutclassdf$psgenes),ylab="Positive selection",color = c("lightblue","red"),main = "Mutation type",cex.axis =1,xlab=paste("p-value(chisq.test)",round(mutclassChisq$p.value,5),sep="="))
mosaicplot(table(gendf$class,gendf$psgenes),ylab="Positive selection",color = c("lightblue","red"),main = "Genetics",cex.axis =1,xlab=paste("p-value(chisq.test)",round(genChisq$p.value,3),sep="="))
mosaicplot(table(mutefdf$class,mutefdf$psgenes),ylab="Positive selection",color = c("lightblue","red"),main = "Mutation type",cex.axis =1,xlab=paste("p-value(chisq.test)",round(mutefChisq$p.value,3),sep="="))
mosaicplot(table(tisdf$class,tisdf$psgenes),ylab="Positive selection",color = c("lightblue","red"),main = "Tissue type",cex.axis =1,xlab=paste("p-value(chisq.test)",round(tisChisq$p.value,3),sep="="))
mosaicplot(table(candf$class,candf$psgenes),ylab="Positive selection",color = c("lightblue","red"),main = "Cancer role",cex.axis =1,xlab=paste("p-value(chisq.test)",round(canChisq$p.value,3),sep="="))
mosaicplot(table(fulldata$driver,fulldata$psgenes),ylab="Positive selection (germinal)",color = c("lightblue","red"),main = "Positive selection (somatic)",cex.axis =1,xlab=paste("p-value(chisq.test)",round(drivChisq$p.value,3),sep="="))

###Correlation dN/dS Somatic vs Germinal evolution
cordNdS<-cor.test(dNdS, log10(dnds_mc))
somgerplot <- ggplot(fulldata,aes(dNdS,log10(pnps),colour=driver))+geom_point(size=2) + xlab("germinal rate (dN/dS)") + ylab("log10(somatic rate (pN/pS))") 
somgerplot <- somgerplot + annotate(x=0.5, y=1.5, label=paste("R = ", round(cordNdS$estimate,2)),geom="text", size=5) + annotate(x=0.5, y=1.35, label=paste("p = ", round(cordNdS$p.value,2)),geom="text", size=5)
somgerplot <- somgerplot + theme(legend.position = c(0.9,0.8), legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=12))
