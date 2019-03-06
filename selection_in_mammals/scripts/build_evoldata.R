library(stringr)
library(magrittr)

genedir<-"/home/uvi/be/avs/cancer_genes_selection/selection_in_mammals/genes"
genelist<-list.files(genedir)
evoldata<-data.frame(gene=genelist)

#Add columns with number of sequences and sequence length
nseqs<-as.numeric()
seqlength<-as.numeric()

for (i in 1:length(genelist)){
  
myfile<-readLines(paste(genedir,"/",genelist[i],"/paml/","M0/out",sep=""))
                  
nseqs[i]<-as.numeric(sapply(strsplit(myfile[4],split = " "),"[",4))
seqlength[i]<-as.numeric(sapply(strsplit(myfile[4],split = " "),"[",8))
}

#Add columns with global dN,dS and dN/dS ratio
dN<-as.numeric()
dS<-as.numeric()
dNdS<-as.numeric()

for (i in 1:length(genelist)){
  
  myfile<-readLines(paste(genedir,"/",genelist[i],"/paml/","M0/out",sep=""))
  
  dN[i]<-as.numeric(str_sub(grep("dN:",myfile,value=T), start=-6))
  dS[i]<-as.numeric(str_sub(grep("dS:",myfile,value=T), start=-7))
  dNdS[i]<-as.numeric(str_sub(grep("omega",myfile,value = T), start=-7))
}


#Add LRTs for test of variation amog sites (M0 vs M3 models in PAML)
lhM0=as.numeric()
lhM3=as.numeric()
for (i in 1:length(genelist)){
  
  myfileM0<-readLines(paste(genedir,"/",genelist[i],"/paml/","M0/out",sep=""))
  lhM0[i]<-as.numeric(sapply(strsplit(myfileM0[grep("lnL",myfileM0)],split=": |    "),"[",4))
  
  myfileM3<-readLines(paste(genedir,"/",genelist[i],"/paml/","M3/out",sep=""))
  lhM3[i]<-as.numeric(sapply(strsplit(myfileM3[grep("lnL",myfileM3)],split=": |    "),"[",4))
  
}

lrtM0M3<-1-pchisq(2*(lhM3-lhM0),6) #Likelihood ratio test

#Add LRTs from the comparisons of PAML site-models
##M1 vs M2 test
lhM1=as.numeric()
lhM2=as.numeric()
p2M2=numeric()
w2M2=numeric()

for (i in 1:length(genelist)){
  
  myfileM1<-readLines(paste(genedir,"/",genelist[i],"/paml/","M1/out",sep=""))
  lhM1[i]<-as.numeric(sapply(strsplit(myfileM1[grep("lnL",myfileM1)],split=": |    "),"[",4))
  
  myfileM2<-readLines(paste(genedir,"/",genelist[i],"/paml/","M2/out",sep=""))
  lhM2[i]<-as.numeric(sapply(strsplit(myfileM2[grep("lnL",myfileM2)],split=": |    "),"[",4))
  
  p2M2[i]<-as.numeric(sapply(strsplit(myfileM2[grep("p:   ",myfileM2)],split="  "),"[",4))
  w2M2[i]<-as.numeric(str_sub(myfileM2[grep("w:  ",myfileM2)], start=-9))
}

p2w2M2<-p2M2*w2M2
lrtM1M2<-1-pchisq(2*(lhM2-lhM1),1) #Likelihood ratio test
padjM1M2<-p.adjust(lrtM1M2,method="fdr", n=length(lrtM1M2)) #Correct p-value for multiple testing


##M8 vs M8a test
lhM8=as.numeric()
lhM8a=as.numeric()


for (i in 1:length(genelist)){
    
  myfileM8a<-readLines(paste(genedir,"/",genelist[i],"/paml/","M8a/out",sep=""))
  lhM8a[i]<-as.numeric(sapply(strsplit(myfileM8a[grep("lnL",myfileM8a)],split=": |    "),"[",4))
  
  myfileM8<-readLines(paste(genedir,"/",genelist[i],"/paml/","M8/out",sep=""))
  lhM8[i]<-as.numeric(sapply(strsplit(myfileM8[grep("lnL",myfileM8)],split=": |    "),"[",4))
  
  
}

lrtM8aM8<-1-pchisq(2*(lhM8-lhM8a),1) #Likelihood ratio test
padjM8aM8<-p.adjust(lrtM8aM8,method="fdr", n=length(lrtM8aM8)) #Correct p-value for multiple testing

#Number of sites under positive selection
pssM2<-as.numeric()
pssM8<-as.numeric()
 
 for (i in 1:length(genelist)){
   myfile_lrtM2<-readLines(paste(genedir,"/",genelist[i],"/paml/",genelist[i],"_test_M1-M2.out",sep=""))
   pssM2[i]<-length(grep("selected",myfile_lrtM2))
   
   myfile_lrtM8<-readLines(paste(genedir,"/",genelist[i],"/paml/",genelist[i],"_test_M7-M8.out",sep=""))
   pssM8[i]<-length(grep("selected",myfile_lrtM8))
 }

#Add the significance (p-value) from LRT of BUSTED analysis
lrtBusted<-numeric()

for (i in 1:length(genelist)) {
  fileBusted<-readLines(paste(genedir,"/",genelist[i],"/hyphy/",genelist[i],"_busted.out",sep=""))
  lrtBusted[i]<-as.numeric(sapply(strsplit(tail(fileBusted,1),split="\\*\\*"),"[",2) %>% 
                             str_sub( .,start=-6))
}

padjBusted<-p.adjust(lrtBusted,method="fdr", n=length(lrtBusted)) #Correct p-value for multiple testing

PS_Busted<-padjBusted< 0.01
PS_Busted<-gsub(TRUE, "Yes",PS_Busted)
PS_Busted<-gsub(FALSE, "No",PS_Busted)

#Conacatenate data
evoldata<-cbind(evoldata,nseqs,seqlength,dN,dS,dNdS,
                lhM1,lhM2,lrtM1M2,padjM1M2, p2w2M2,
                lhM8a,lhM8,lrtM8aM8, padjM8aM8, pssM2, pssM8,
                lrtBusted,padjBusted)
