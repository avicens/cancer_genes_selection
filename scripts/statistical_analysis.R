#Combine functional and evolutionary data frames
setwd("Dropbox/phylogenomics_lab_dbx/cancer_genes_selection/") #Set working directory
load("data/evoldata") #Load data frame with evolutionary data (PAML)
load("data/fundata") #Load data frame with functional data (COSMIC)

fulldata<-merge(evoldata, fundata,by="gene")
attach(fulldata)

### Comparison of dN/dS across categories
source("scripts/dndspercat.R") #Load function to generate dataframes per category
source("scripts/make_comparisons.R") #Load function to generate comparisons invariables with more than 2 classes
library(ggpubr) #Load "ggpubr" library for getting labelled box plot

###Statistical analysis

#1.Mutations type (somatic, germline)
mutclasses<-unique(unlist(strsplit(as.character(fulldata$Mut.class),split=", ")))
mutclassnames<-c("Somatic","Som + Germ", "Germinal")
mutclassdf<-dndspercat(mutclasses,mutclassnames,fulldata,Mut.class)

mutclasscomp<-make.comparisons(mutclassnames) #Generate a list with all possible comparisons

#Boxplot
mut_dNdS<-ggboxplot(mutclassdf, x = "class", y = "log10(dNdS)", select= mutclassnames, fill="class",palette="jco",xlab=FALSE, ylab = "log10 (dN/dS)") 
mut_dNdS2 <- mut_dNdS + stat_compare_means(method = "t.test", comparisons = mutclasscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 1.5, label.x = 0.7) + theme(legend.position="none")

##2. Molecular genetics
genet<-levels(fulldata$inheritance)
genet<-c("Dom","Rec")
genames<-c("Dominant","Recessive")
gendf<-dndspercat(genet,genet,fulldata,inheritance)

#Boxplot
gencomp<-make.comparisons(genames) #Generate a list with all possible comparisons

inh_dnds<-ggboxplot(fulldata,x="inheritance",y="log10(dNdS)", fill ="inheritance", palette="jco", select=c("Dom","Rec"),xlab=FALSE, ylab="log10(dN/dS)")
inh_dnds2 <-inh_dnds + stat_compare_means(method = "t.test", label.y = 1, label.x = 0.7, cex=5) + theme(legend.position="none")

#3. Functional effect of coding mutations
mutype<-gsub("\"","",fulldata$impact)
mutype<-unique(unlist(strsplit(mutype,", |. |,| ")))
mutype2<-c("Mis","N")
mutnames<-c("Missense","Nonsense")
mutefdf<-dndspercat(mutype2,mutnames,fulldata,impact)

#Boxplot
imp_dnds<-ggboxplot(mutefdf, x = "class", y = "log10(dNdS)", select= mutnames, fill="class",palette="jco",xlab=FALSE, ylab = "Global dN/dS") 
mutef_dnds2<- mutef_dnds + stat_compare_means(method = "t.test", label.y = 0.6, label.x = 1.5, cex= 5)

##4. Tissue type
tissues<-unique(unlist(strsplit(levels(fulldata$tissue.type),split=",|;")))
tisnames<-c("Epith","Leu/Lym","Mesen","Other")
tisdf<-dndspercat(tissues,tisnames,fulldata,tissue.type)

#Data exploration
table(tisdf$class) #Check the sample size for each group
tiscomp<-make.comparisons(tisnames) #Generate a list with all possible comparisons

#Boxplot
#dN/dS
tis_dnds<-ggboxplot(tisdf, x = "class", y = "log10(dNdS)", select= tisnames, fill="class",palette="jco",xlab=FALSE, ylab = "log10 (dN/dS)") 
tis_dnds2<-tis_dnds + stat_compare_means(method = "t.test", comparisons = tiscomp,label="p.signif") + stat_compare_means(method = "anova", label.y = 1.2, label.x = 0.7) + theme(legend.position="none")

##5.Cancer role
canrole<-gsub("\"","",levels(fulldata$cancer.role))
canrole<-unique(unlist(strsplit(canrole,", ")))
candf<-dndspercat(canrole, canrole,fulldata,cancer.role)

cancomp<-make.comparisons(canrole) #Generate a list with all possible comparisons

#dN/dS
can_dnds<-ggboxplot(candf, x = "class", y = "log10(dNdS)", select= canrole, fill="class",palette="jco",xlab=FALSE, ylab = "log10 (dN/dS)") 
can_dnds2 <- can_dnds + stat_compare_means(method = "t.test", comparisons = cancomp, label ="p.signif") + stat_compare_means(method = "anova", label.y = 1, label.x = 0.7, cex=5) + theme(legend.position = "none")

##6.Chromosome (X/autosomic)
chr_dNdS<-ggboxplot(fulldata,x="chr.type",y="log10(dNdS)", fill="chr.type",palette="jco", xlab=FALSE,ylab="log10(dN/dS)")
chr_dNdS2 <- chr_dNdS + stat_compare_means(method = "t.test", label.y = 1, label.x = 0.7) + theme(legend.position = "none")

##multiplot (dN/dS)
ggarrange(mutclass_dNdS2, gendf, mutef_dnds2, tis_dnds2, can_dnds2, driv_dNdS2, nrow=2, ncol=3)

##Multiplot (dN,dS)
ggarrange(mutclass_dN2,mutclass_dS2,gen_dN2,gen_dS2,mutef_dN2, mutef_dS2,tis_dN2,tis_dS2,can_dN2,can_dS2, driv_dN2, driv_dS2, nrow=2, ncol= 6)

#Comparison of proportion of positively selected genes (chi-square text)
mutclassChisq<-chisq.test(table(mutclassdf$PS, mutclassdf$class))
genChisq<-chisq.test(table(gendf$psgenes, gendf$class))
mutefChisq<-chisq.test(table(mutefdf$psgenes, mutefdf$class))
tisChisq<-chisq.test(table(tisdf$psgenes, tisdf$class))
canChisq<-chisq.test(table(candf$psgenes, candf$class))
drivChisq<-chisq.test(table(fulldata$psgenes,fulldata$driver))

mut_ps<-as.data.frame(table(mutclassdf$PS,mutclassdf$class))
names(mut_ps)<-c("PS","mut.class","Freq")
mut_ps2<-ggplot(mut_ps,aes(mut.class,Freq,fill=PS)) + geom_bar(stat="identity",position = position_fill(reverse = TRUE)) + labs(x="Mutation Type", y="Proportion of genes") + scale_fill_manual(values = c("DodgerBlue","FireBrick")) + theme_grey(base_size=14) + scale_x_discrete(limits=c("Somatic","Som + Germ", "Germinal")) + theme(legend.position ="top")

inh_ps<-as.data.frame(table(PS,inheritance))
inh_ps2<-ggplot(inh_ps,aes(inheritance,Freq,fill=PS)) + geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + labs(x="Inheritance", y="Proportion of genes") + scale_fill_manual(values = c("DodgerBlue","FireBrick")) + theme_grey(base_size=14) + scale_x_discrete(limits=c("Dom","Rec")) + theme(legend.position ="top")

tis_ps<-as.data.frame(table(tisdf$PS,tisdf$class))
names(tis_ps)<-c("PS","tissue.type","Freq")
tis_ps2<-ggplot(tis_ps,aes(tissue.type,Freq,fill=PS)) + geom_bar(stat="identity",position = position_fill(reverse = TRUE)) + labs(x="Tissue Type", y="Proportion of genes") + scale_fill_manual(values = c("DodgerBlue","FireBrick")) + theme_grey(base_size=14) + theme(legend.position ="top")

can_ps<-as.data.frame(table(candf$PS,candf$class))
names(can_ps)<-c("PS","cancer.role","Freq")
can_ps2<-ggplot(can_ps,aes(cancer.role,Freq,fill=PS)) + geom_bar(stat="identity",position = position_fill(reverse = TRUE)) + labs(x="Role in Cancer", y="Proportion of genes") + scale_fill_manual(values = c("DodgerBlue","FireBrick")) + theme_grey(base_size=14) + theme(legend.position ="top")

chr_ps<-as.data.frame(table(PS,chr.type))
chr_ps2<-ggplot(chr_ps,aes(chr.type,Freq,fill=PS)) + geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + labs(x="Chromosome", y="Number of genes") + scale_fill_manual(values = c("DodgerBlue","FireBrick")) + theme_grey(base_size=14) + scale_x_discrete(limits=c("A","X")) + theme(legend.position ="top")

png("statistical_comparison/Rplot02.png",width = 1300, height=650)
ggarrange(mut_dNdS2,inh_dnds2, tis_dnds2, can_dnds2, chr_dNdS2,mut_ps2,inh_ps2, tis_ps2, can_ps2, chr_ps2,ncol = 5, nrow=2)
dev.off()

###Correlation dN/dS Somatic vs Germinal evolution
cordNdS<-cor.test(dNdS, log10(dnds_mc))
somgerplot <- ggplot(fulldata,aes(dNdS,log10(pnps),colour=driver))+geom_point(size=2) + xlab("germinal rate (dN/dS)") + ylab("log10(somatic rate (pN/pS))") 
somgerplot <- somgerplot + annotate(x=0.5, y=1.5, label=paste("R = ", round(cordNdS$estimate,2)),geom="text", size=5) + annotate(x=0.5, y=1.35, label=paste("p = ", round(cordNdS$p.value,2)),geom="text", size=5)
somgerplot <- somgerplot + theme(legend.position = c(0.9,0.8), legend.title = element_text(size=14, face="bold"), legend.text = element_text(size=12))
