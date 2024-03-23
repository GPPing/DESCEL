library(poolr)  
library(metap)
library(SingleCellExperiment)
library(BPSC)
library(DEsingle)
library(DESeq2)
library(MAST)
library(monocle)
library(scDD)
library(limma)
library(zingeR)

#load example data
#input gene expression matrix and cell label information
setwd("C:/Users/ping/Desktop/DESCEL-main/data")
load("count.RData")
load("group.RData")

#set your working path
setwd("C:/Users/ping/Desktop/DESCEL-main/R")
source("DESCEL.R")

#data pre-processing and differential expression analysis
#Here we can choose different normalization methods, such as CPM, TMM, RLE, the default setting is CPM



Pvals <- DESCEL(raw.count = count, cell.label = group)
p.ACAT <- p.ACAT(Pvals, weights = TRUE)
#p.ACAT <- combination.Pvals(Pvals, weights = FALSE)
FDR <- DESCEL.p.adjust(p.ACAT, adjusted.method = "fdr")




