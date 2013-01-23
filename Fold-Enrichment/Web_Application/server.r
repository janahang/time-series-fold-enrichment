# runApp("~/Web_Application/")
rm(list=ls())
library("fmsb")
library(ggplot2)
library(shiny)
source("~\\Fold_Enrichment\\Web_Application\\radar_chart.r")
setwd("~\\Fold_Enrichment")

# Since they used two different microarrays in their experiment, 
# Common probesets were chosen for calculationg fold enrichment
u133a = read.table("GPL3921_expression.txt",sep="\t",header=TRUE,row.names=1)
u133p2 = read.table("GPL570_expression.txt",sep="\t",header=TRUE,row.names=1) 
pbmc_subset = read.table("PBMC_subset_data.txt",sep="\t",header=TRUE,row.names=1)
bcell_subset = read.table("B_cell_subset_data.txt",sep="\t",header=TRUE,row.names=1)

sample = read.table("Sample_input.txt",sep="\t",header=TRUE,row.names=1)
gene_list = sample
gene_list = gene_list[complete.cases(gene_list[,1]),]
gene_list = gene_list[!duplicated(gene_list$gene_symbol),]
rownames(gene_list) = gene_list$gene_symbol

#######################################################
# Obatain the list of genes for each subset of interest
#######################################################
# PBMC Subsets
#######################################################
bcell = pbmc_subset[grep("B-cell",pbmc_subset$Highly.Expressed.in.PBMC.type),]
bcell = bcell[complete.cases(bcell[,1]),]
bcell = bcell[-(which(bcell$gene_symbol=="abParts")),]
bcell = bcell[!duplicated(bcell$gene_symbol),]
bcell = bcell[,1]

tcell = pbmc_subset[grep("T-cell",pbmc_subset$Highly.Expressed.in.PBMC.type),]
tcell = tcell[complete.cases(tcell[,1]),]
tcell = tcell[!duplicated(tcell$gene_symbol),]
tcell = tcell[,1]

mDC = pbmc_subset[grep("mDC",pbmc_subset$Highly.Expressed.in.PBMC.type),]
mDC = mDC[complete.cases(mDC[,1]),]
mDC = mDC[!duplicated(mDC$gene_symbol),]
mDC = mDC[,1]

pDC = pbmc_subset[grep("pDC",pbmc_subset$Highly.Expressed.in.PBMC.type),]
pDC = pDC[complete.cases(pDC[,1]),]
pDC = pDC[-(which(pDC$gene_symbol=="abParts")),]
pDC = pDC[!duplicated(pDC$gene_symbol),]
pDC = pDC[,1]


monocytes = pbmc_subset[grep("Monocytes",pbmc_subset$Highly.Expressed.in.PBMC.type),]
monocytes = monocytes[complete.cases(monocytes[,1]),]
monocytes = monocytes[!duplicated(monocytes$gene_symbol),]
monocytes = monocytes[,1]

NK = pbmc_subset[grep("NK",pbmc_subset$Highly.Expressed.in.PBMC.type),]
NK = NK[complete.cases(NK[,1]),]
NK = NK[-(which(NK$gene_symbol=="abParts")),]
NK = NK[!duplicated(NK$gene_symbol),]
NK = NK[,1]

#######################################################
# B-Cell Subsets
#######################################################
Plasma = bcell_subset[grep("Plasma",bcell_subset$Highly.Expressed.in.B.cells..PBMC..and.highly.expressed.in.B.cell.subsets),]
Plasma = Plasma[complete.cases(Plasma[,1]),]
Plasma = Plasma[-(which(Plasma$gene_symbol=="abParts")),]
Plasma = Plasma[!duplicated(Plasma$gene_symbol),]
Plasma = Plasma[,1]


Memory = bcell_subset[grep("Memory",bcell_subset$Highly.Expressed.in.B.cells..PBMC..and.highly.expressed.in.B.cell.subsets),]
Memory = Memory[complete.cases(Memory[,1]),]
Memory = Memory[!duplicated(Memory$gene_symbol),]
Memory = Memory[,1]

Naive = bcell_subset[grep("Naive",bcell_subset$Highly.Expressed.in.B.cells..PBMC..and.highly.expressed.in.B.cell.subsets),]
Naive = Naive[complete.cases(Naive[,1]),]
Naive = Naive[!duplicated(Naive$gene_symbol),]
Naive = Naive[,1]


GC = bcell_subset[grep("GC",bcell_subset$Highly.Expressed.in.B.cells..PBMC..and.highly.expressed.in.B.cell.subsets),]
GC = GC[complete.cases(GC[,1]),]
GC = GC[-(which(GC$gene_symbol=="abParts")),]
GC = GC[!duplicated(GC$gene_symbol),]
GC = GC[,1]

gene_list = sample
pbmc = list(bcell,tcell,mDC,pDC,monocytes,NK)
results_p = NULL
results_p = matrix(nrow=ncol(gene_list),ncol=6)
rownames(results_p) = colnames(gene_list)
colnames(results_p) = c("B-cell","T-cell","mDC","pDC","Monocytes","NK cell")

bcell_subset = list(Plasma,Memory,Naive,GC)
results_b = NULL
results_b = matrix(nrow=ncol(gene_list),ncol=4)
rownames(results_b) = colnames(gene_list)
colnames(results_b) = c("Plasma","Memory","Naive","GC")


# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
  
  # Function that generates a plot of the distribution. The function
  # is wrapped in a call to reactivePlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically 
  #     re-executed when inputs change
  #  2) Its output type is a plot 
  #
  output$distPlot <- reactivePlot(function() {
    
    
    for(j in 1:ncol(gene_list)){
      gene_list = gene_list[which(gene_list[,j]==input$condition),]
      for(i in 1:length(pbmc)){
        zx = length(intersect(rownames(gene_list),pbmc[[i]]))
        deg_z = length(rownames(gene_list))
        subset_x = length(pbmc[[i]])
        total = length(intersect(rownames(u133p2),rownames(u133a)))
        fold_enrichment = (zx/deg_z)/(subset_x/total)
        results_p[j,i] = fold_enrichment
      }
      for(i in 1:length(bcell_subset)){
        zx = length(intersect(rownames(gene_list),bcell_subset[[i]]))
        deg_z = length(rownames(gene_list))
        subset_x = length(bcell_subset[[i]])
        total = length(intersect(rownames(u133p2),rownames(u133a)))
        fold_enrichment = (zx/deg_z)/(subset_x/total)
        results_b[j,i] = fold_enrichment
      }
    }    
    
    results_p[is.nan(results_p)] <- 0 
    results_b[is.nan(results_b)] <- 0 
    par(mfrow=c(2,1))
    max = ceiling(max(results_p)*1.5)
    min = floor(min(results_p)*0.5)
    x = c(rep(max,6),rep(min,6))
    maxmin = as.data.frame(matrix(x,nrow=2, ncol=6,byrow=T))
    colnames(maxmin) = c("B-cell","T-cell","mDC","pDC","Monocytes","NK cell")
    data_p = as.data.frame(rbind(maxmin,results_p))
    datap = data_p[-(1:2),]
    color = c(2:(nrow(datap)+1))
    radarchart(data_p,axistype=0,plty=1,plwd=3,pcol=color,maxmin=TRUE,
               cglty=1,cglcol=1,
               labels = seq(from = min(x), to = max(x), length = 5),
               title="Gene Enrichment for PBMC subtypes")
    legend(-3,1.1,rownames(datap),col=color, pch=20,cex=0.8,border="white",lty=2, lwd=4,bty='n')
    max = ceiling(max(results_b)*1.5)
    min = floor(min(results_b)*0.5)
    x = c(rep(max,4),rep(min,4))
    maxmin = as.data.frame(matrix(x,nrow=2, ncol=4,byrow=T))
    colnames(maxmin) = c("Plasma","Memory","Naive","GC")
    data_b = as.data.frame(rbind(maxmin,results_b))
    datab = data_b[-(1:2),]
    seg = abs(min)+abs(max)+1
    radarchart(data_b,axistype=0,plty=1,plwd=3,pcol=color,maxmin=TRUE,
               cglty=1,cglcol=1,
               labels = seq(from = min(x), to = max(x), length = 5),
               title="Gene Enrichment B-cell subtypes")
     })
})
