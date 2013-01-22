Fold-Enrichment
===============

Create a Shiny web application using R for Fold Enrichment Analyses of a Time-Series Experiment

Fold Enrichment Analyses for a Time-Series Experiment

This project involves reproducing the analysis and spider plot from Figure 3b,c,e in (Nakaya H et al., 2011,http://www.ncbi.nlm.nih.gov/pubmed/21743478). 
These plots show the enrichment of cell-specific gene signatures in their datasets. 
The project would be to write an R Shiny web application that recreates their analysis so it can be applied to any other data sets. 
The supplementary methods in the paper go into some detail about exactly how they calculate the "fold enrichment" statistic. 
The cell specific gene sets are available in their supplementary materials. 
The input to the function would be a file with a list of gene symbols (rows) with differential expression condition (up-regulated,down-regulated,no change) under each time of experiment (columns). 
Each time-series would appear as a separate line on the spider plot.

###############################
# Instructions
###############################
1) The input file should have the columns "genes","time series 1","time series 3","time series 3".. ,"time series n"
2) The differential expression conditions should represent as 1=up-regulated, 0=unchanged, -1=down-regulated.
3) The code could be modified according to your needs.
4) Set your working directory properly
5) In order to run the app, from your R teminal type 'runApp("~/Web_Application/")' [change the path accordingly]
