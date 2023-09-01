library(ggpubr)
library(DescTools)
library(vegan)
library(immunarch)
library(tidyverse)
library(data.table)

getdiv <- function(immdata) {
    TRB <- data.frame()
    for (n in immdata$meta$Sample){
        shannon_entropy <- vegan::diversity(immdata$data[[n]]$Clones , "shannon")
        Pielou_index <- shannon_entropy / log(vegan::specnumber(immdata$data[[n]]$Clones))
        Clonality = 1 - Pielou_index
        Gini.Simpson <- vegan::diversity(immdata$data[[n]]$Clones, "simpson")
        inverse.Simpson <- vegan::diversity(immdata$data[[n]]$Clones, "invsimpson")

        d50 <- 0
        i <- 1
        lst <- c()
        while (d50 < 0.5) {
            lst <- append(lst,i)
            d50 <- d50 + immdata$data[[n]][i,"Proportion"]
            i=i+1
        }

        d50 <- lst[length(lst)] 
        DE50 <- d50/nrow(immdata$data[[n]])

        diversity <- data.frame(
            sample = n,
            unique.clonotypes = nrow(immdata$data[[n]]),
            total.clonotypes = sum(immdata$data[[n]]$Clones),
            shannon = shannon_entropy,
            Gini_Simpson = Gini.Simpson,
            inv_Simpson = inverse.Simpson,
            Pielou_index = Pielou_index,
            Clonality = Clonality,
            d50 = d50,
            DE50 = DE50)
            
    TRB <- rbind(TRB, diversity)
}
  return(TRB)
}
