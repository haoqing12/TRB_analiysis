# 1. The consensus of local nucleotide and amino acid sequences in the CDR3 region was analyzed by aligning and calculating sequence similarity with the “msa” package.
library(msa)
mySequenceFile <- system.file("examples", "aa.fasta", package="msa") 
mySequences <- readAAStringSet(mySequenceFile) 
mySequences
myFirstAlignment <- msa(mySequences)
myFirstAlignment
print(myFirstAlignment, show="complete")
msaPrettyPrint(myFirstAlignment, output="asis",  showNames="none", showLogo="top", logoColors="rasmol", shadingMode="similar", showConsensus= "bottom", consensusColors= "GreenRed",showLegend=FALSE, askForOverwrite=FALSE)

# 2. The consensus of local nucleotide and amino acid sequences in the CDR3 region was plotted with the “ggseqlogo” packages.
library(ggplot2)
library(ggseqlogo)
rm(list=ls())
setwd('C:\\Users\\paper\\Desktop\\TCRseq')
motif_seq=read.delim(header=F,file='logo_demo.txt',sep='\t',stringsAsFactors = F,na.strings = '-')
g=ggseqlogo(motif_seq$V1, method="bits")
g

# 3. Differential analysis was performed using 'limma' package.

library("limma")
setwd("C:\\Users\\paper\\Desktop\\TCRseq") 
inputFile="symbol.txt" 
fdrFilter=0.05   
logFCfilter=0.263  
conNum=19                                                       
treatNum=19  
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
write.table(outTab,file="all.txt",sep="\t",row.names=F,quote=F)

# 4. Volcano plot
library(ggplot2)
setwd('C:\\Users\\paper\\Desktop\\TCRseq')
temp = read.csv("all.csv",header=T,row.names=1)
p=ggplot(temp,aes(logFC,-1*log10(Pvalue)))
p+geom_point()
p + geom_point(color ="red")
p +geom_point(aes(color ="red"))
p + geom_point(aes(color =Significant))
r03xy = p +geom_point(aes(color =Significant)) + xlim(-6,6) + ylim(0,4)
r03xy + labs(title="Volcanoplot",x="log2(FC)")
r03xyp = r03xy + labs(title="Volcanoplot",x="log2(FC)")
r03xyp + scale_color_manual(values =c("royalblue3","black", "firebrick3"))
volcano=r03xyp + scale_color_manual(values =c("royalblue3","black", "firebrick3"))
volcano+geom_hline(yintercept=1.3)+geom_vline(xintercept=c(-0.263,0.263))
volcano+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-0.263,0.263),linetype=4)

# 5. Heatmap
library(pheatmap)
setwd("C:\\Users\\paper\\Desktop\\TCRseq")         
data<- read.table("heatmap.txt",head=T,row.names = 1,sep="\t")
pheatmap(data,cluster_cols = F,cluster_rows = F,show_colnames=F,show_rownames = T,border_color = NA)

# 6. Paired box plot
library(ggplot2)
library(ggpubr)
setwd("C:\\Users\\paper\\Desktop\\TCRseq") 
cv=read.table("box.txt",header=T,sep="\t")
cvbox=ggplot(cv,aes(x=TCR,y=Frequency,fill=treatment))+geom_boxplot(outlier.size = -1)
cvbox+geom_jitter() 
compare_means(Frequency~treatment,data=cv,method = "wilcox.test")
my_comparison=list(c("A","B"))
cvbox+stat_compare_means(comparisons = my_comparison)+ 
  stat_compare_means(method = "wilcox.test") +geom_jitter() 

