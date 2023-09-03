limma.R 
##Differential gene expression analysis
library(limma)
setwd("C:\\Users\\Administrator\\Desktop\\project")
#normalize
rt=read.table("mRNA_FPKM.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt=normalizeBetweenArrays(as.matrix(rt))
rt=log2(rt+1)

Type=c(rep("con",16),rep("treat",21))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

all=topTable(fit2,number=Inf,adjust.method="holm")
all=rbind(id=colnames(all),all)
write.table(all,file="all.txt",sep="\t",quote=F,col.names=F)
diffSig = all[(all$P.Value < PValue & (all$logFC>logFCcutoff | all$logFC<(-logFCcutoff))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = all[(all$P.Value < PValue & (all$logFC>logFCcutoff)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = all[(all$P.Value < PValue & (all$logFC<(-logFCcutoff))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

#write expression level of diff gene
diffName=row.names(diffSig)
hmExp=rt[diffName,]
hmExp=rbind(id=colnames(hmExp),hmExp)
write.table(hmExp,file="heatmap.txt",sep="\t",quote=F,col.names=F)

#volcano
library(ggplot2)
library(ggrepel)
data = read.table("all.txt", header=TRUE)

####logFC为1.5
m=ggplot(data=data, aes(x=logFC, y =-log10(P.Value))) +
  geom_point(data=subset(data,abs(data$logFC) <= 1.5),aes(size=abs(logFC)),color="gray",alpha=0.5) +
  geom_point(data=subset(data,data$P.Value >= 0.05&abs(data$logFC) > 1.5),aes(size=abs(logFC)),color="gray",alpha=0.5) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC > 1.5),aes(size=abs(logFC)),color="red",alpha=0.5) +
  geom_point(data=subset(data,data$P.Value<0.05 & data$logFC < -1.5),aes(size=abs(logFC)),color="darkgreen",alpha=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  labs(x="Log2FC",y="-Log10 (P Value)")+
  theme(legend.position='none')
#geom_text_repel(data=subset(data, abs(logFC) > 5), aes(label=id),color="black",alpha = 0.8)
print(m) 
ggsave(m,file="volcano.pdf",width = 8,height = 6)

#Note: Advanced volcano plot was performed using the OmicStudio tools at https://www.omicstudio.cn/tool.

#Heatmap
library(pheatmap)
setwd("C:\\Users\\paper\\Desktop\\TCRseq")         
data<- read.table("heatmap.txt",head=T,row.names = 1,sep="\t")
pheatmap(data,cluster_cols = F,cluster_rows = F,show_colnames=F,show_rownames = T,border_color = NA)


#GO enrichment analysis

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

setwd("C:\\Users\\Administrator\\Desktop\\project")

library("org.Hs.eg.db") 
rt=read.table("diff.txt",sep="\t",check.names=F,header=T) 
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="entrezID.txt",sep="\t",quote=F,row.names=F)

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("enrichplot")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("C:\\Users\\Administrator\\Desktop\\project")
rt=read.table("entrezID.txt",sep="\t",header=T,check.names=F)   
rt=rt[is.na(rt[,"entrezID"])==F,]                           
gene=rt$entrezID

kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="BP",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)
dotplot(kk,showCategory = 10) 

pdf(file="barplot.pdf",width = 10,height = 8)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#bubbleplot
library(ggplot2)
GO = read.csv("C:\\Users\\Administrator\\Desktop\\project\\down2.csv",header=TRUE,row.names=1,check.names = FALSE)   
p = ggplot(GO,aes(Pvalue,GO))
p=p + geom_point()  
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=-1*log10(Pvalue)))
pr = pbubble+scale_color_gradient(low="gold2",high = "royalblue3")
pr = pr+labs(color=expression(-log[10](Pvalue)),size="Count",  
             x="P value",y="Biological process",title="Biological process")
pr 


##ESTIMATE
##library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

library(estimate)
setwd("C:\\Users\\Desktop\\ESTIMATE")

filterCommonGenes(input.f="symbol_T365.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="ESTIMATE_scores.txt",sep="\t",quote=F,col.names=F)

 CIBERSORT R script v1.03
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.pacakges('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
#' single-threaded in Windows.
#'
#' Usage:
#'       Navigate to directory containing R script
#'
#'   In R:
#'       source('CIBERSORT.R')
#'       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#'
#'       Options:
#'       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#'       ii) QN = Quantile normalization of input mixture (default = TRUE)
#'
#' Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
#' Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' Core algorithm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
CoreAlg <- function(X, y){

  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#' do permutations
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#' Main functions
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mixture_file heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @export
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  library(e1071)
  library(parallel)
  library(preprocessCore)

  #read in data
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }

  #save results
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "symbol_T365.txt", perm=100, QN=FALSE)

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(pheatmap)
library(ggcorrplot)

# Read TME file
TME.results <-read.table("CIBERSORT-Results.txt",sep = "\t",header = T,row.names = 1)

# Read phenotype file of TCGA
pheno.data <- read.table("phenotype.txt", header = T, sep = "\t")

# Merge TME file with phenotype data
sample.list <- as.character(unique(pheno.data$PATIENT_ID))
pheno.data <-
  pheno.data[match(sample.list, pheno.data$PATIENT_ID), ]
TME.results <-
  TME.results[match(sample.list, rownames(TME.results)), ]
TME.data <- cbind(pheno.data, TME.results)

# Show TME cells
TME.cells <- colnames(TME.results)[1:22]
TME.cells

#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(
    PATIENT_ID = TME.data$PATIENT_ID,
    CellType = TME.cells[i],
    response = TME.data$response,
    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}

#绠卞紡鍥綽oxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#鏉″舰鍥綛arplot of cell componment of each sample using hclust
sample.index <-
  hclust(dist(TME.data[, TME.cells]), method = "ward.D")$order
sample.order <- TME.data$PATIENT_ID[sample.index]

ggbarplot(
  plot.info,
  x = "PATIENT_ID",
  y = "Composition",
  size = 0,
  fill = "CellType",
  color = "CellType",
  order = sample.order
) +
  theme_base() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 1
    ),
    legend.position = "bottom"
  )

# boxplot group by risk
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "response",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  stat_compare_means(
    label = "p.signif",
    method = "t.test",
    ref.group = ".all.",# Pairwise comparison against all
    hide.ns = F
    )+
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))


