
library(stringr)
library(dplyr)
wd<-getwd()


build.node.table<-function(sample) {
    
    id.file<-readLines(paste(wd,"/data/paml/",sample,"/",sample,"_node_ids.txt",sep=""))

    #Build node Ids table
    header.lines<-grep("Node ID :",id.file)
    id.lines<-id.file[(header.lines[1] +2) : (header.lines[2] - 2)]
    s<-as.integer(sapply(strsplit(id.lines,":"),"[",1))
    t<-sapply(strsplit(id.lines,":"),"[",2)
    sp<-paste(sapply(strsplit(t,"_"),"[",c(1)),sapply(strsplit(t,"_"),"[",c(2)))
    id<-sapply(strsplit(t,"_"),"[",3)

    id.df<-data.frame(Node = s, Species = sp, ID = id, stringsAsFactors=F)

    return(id.df)
}

build.lrt.table<- function(sample) {
    fb.file<-readLines(paste(wd,"/data/paml/",sample,"/",sample,"_fb_primates_leaves.out",sep=""))
    
    #Parsing table LRT
    lrt.m0fb<-grep("M0~[a-z0-9]*",fb.file, value=T); lrt.m0fb<-lrt.m0fb[3:length(lrt.m0fb)]
    lh.fb<-sapply(strsplit(lrt.m0fb,"\\|"),"[",2)
    nodes<-as.integer(sapply(strsplit(lh.fb,"\\.|-"),"[",2))
    p.val<-sapply(strsplit(lrt.m0fb,"\\|"),"[",3); p.val<-as.numeric(gsub(" |\\*","",p.val))
    p.adj<-p.adjust(p.val, n=length(p.val),method="fdr") 
    lrt.df<-data.frame(Node = nodes, LRT = p.val, P.adj = p.adj)


    return(lrt.df)
}

build.rates.table <- function(sample) {

    fb.file<-readLines(paste(wd,"/data/paml/",sample,"/",sample,"_fb_primates_leaves.out",sep=""))

    nodes<-integer()
    bg.omega<-numeric()
    fg.omega<-numeric()

    #Parsing table dN/dS
    ids<-grep("Model b_free.[0-9]+-1$",fb.file)
    for (id in ids) {
        nodes<-append(nodes,as.integer(gsub(" - Model b_free.|-1", "", fb.file[id])))
        omegas<-fb.file[c(id+5, id+6)]

        bg.omega<-append(bg.omega,as.numeric(sapply(strsplit(grep("background",omegas, value =T),"=>"),"[",2)))
        fg.omega<-append(fg.omega,as.numeric(sapply(strsplit(grep("#1",omegas, value =T),"=>"),"[",2)))
    }

    omega.df<-data.frame(Node = nodes, Bg.omega = bg.omega, Fg.omega = fg.omega)

}

samples.list<-list.files("data/paml/", recursive=T, pattern="fb_primates_leaves.out")
#samples<-sapply(strsplit(samples.list,"/"),"[",1)
samples<-read.table("genelists/cancer_genes_signif_human_dnds_EOGid.txt",header=F,stringsAsFactors=F)[,1]

fb.df<-data.frame(Node=integer(), Species=character(), Id=character(), LRT = numeric(), P.adj = numeric(),Bg.omega = numeric(), Fg.omega = numeric(),Group=character())

for (smp in samples) {
    cat("Attaching data for group",smp,"\n")
    
    smp.nodes<-build.node.table(smp)
    smp.lrt<-build.lrt.table(smp)
    smp.rates<-build.rates.table(smp)

    smp.df<-Reduce(function(x,y) merge(x,y),list(smp.nodes, smp.lrt, smp.rates))

    smp.df$Group<-rep(smp,nrow(smp.df))

    fb.df<-rbind(smp.df)

}


