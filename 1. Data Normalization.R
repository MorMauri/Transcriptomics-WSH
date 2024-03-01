#####################################################################################################################################
#####################################################################################################################################
#################################################                                   #################################################
#################################################     N O R M A L I Z A T I O N     ################################################# 
#################################################                                   #################################################
#####################################################################################################################################
#####################################################################################################################################


library(limma)
library(edgeR)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(knitr)
library(piano)
library(pheatmap)
library(org.Sc.sgd.db)
library(shiny)
library(plotly)
library(EnhancedVolcano)
library(grid)
library(tibble)
library(dplyr)
library(data.table)
library(writexl)
library(devtools)
library(ggpubr)
library(eulerr)
library(snowfall)


count.data <- read.delim('DEanalysisMauri/gene_counts/gene_counts_anaerobic.txt', stringsAsFactors=F, row.names=1)

# preview the first few rows of the data
head(count.data)

#you may have a look at the count data in Rstudio
View(count.data)    
nrow(count.data)



#We want to determine which genes are differentially expressed between two groups of samples.
#in this case, between “pheromone-unaffected” and “pheromone-affected” samples, which we have indicated by column names beginning with a “U” or “P”, respectively.
#We can take a quick look at the column names to see what samples we have
colnames(count.data)

#Generate a factor specifying the group to which each sample belongs
sample.groups <- factor(startsWith(colnames(count.data), 'CEN'), labels=c('strain','control'))

sample.groups <- factor(colnames(count.data), levels = c('CEN_ana_1','CEN_ana_2','CEN_ana_3','CEN_ana_4','KE_ana_1','KE_ana_2','KE_ana_3','KE_ana_4','EtOH_ana_1','EtOH_ana_2','EtOH_ana_3','EtOH_ana_4','WT31_ana_1','WT31_ana_2','WT31_ana_3','WT31_ana_4','WT109_ana_1','WT109_ana_2','WT109_ana_3','WT109_ana_4'), 
                        labels=c('CEN.PK','CEN.PK','CEN.PK','CEN.PK','KE6','KE6','KE6','KE6','EtOH','EtOH','EtOH','EtOH','WT31','WT31','WT31','WT31','WT109','WT109','WT109','WT109'))

#check that we set up the groups correctly
data.frame(sample=colnames(count.data), groups=sample.groups)

#Create a DGEList object from the count.data dataframe and sample.groups factor:
y <- DGEList(counts=count.data, group=sample.groups)

#This DGEList object ‘holds’ the dataset to be analysed by edgeR and the subsequent calculations relevant for the analysis.
#Before performing any DE analysis, we can take a look at similarities and differences among different samples using a multidimensional scaling (MDS) plot:
MDS <- plotMDS(y, cex=1.8) #"ces" sets the label font size
plot(MDS, col = factor(sample.groups))
legend("topleft",
       legend = levels(factor(sample.groups)),
       pch = 19,
       col = factor(levels(factor(sample.groups))))
#Often, RNA Seq data will contain low-count genes that are not expressed, or expressed at very low levels in many samples. The data for these genes is essentially noise that does not provide any meaningful information and should therefore be removed prior to DE analysis.
#Although there is no single “correct” approach for removing low-count genes, a common technique involves removing those that are below some threshold in a given fraction of the samples.
#Filter (remove) low-count genes from the count data, using an arbitrary cut off of 10 reads:
#Normalize for library size, by taking the log count-per-million. It's not a satisfactory normalization, but we just use it to plot raw reads
cpm <- cpm(y)
lcpm <- cpm(y, log = T)

#Plot raw logCPM reads before filtering
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired") # Ignore warning, some samples will have identical color.
par(mfrow = c(1, 2)) #This tells to plot charts in 1 line with 2 columns (basically this chart and next are gonna be plotted in the same fig, one on the left one on the right)
#Select this all loop and run it. if you run it line by line it won't work
plot(
  density(lcpm[, 1]),
  col = col[1],
  ylim = c(0, 0.6),
  las = 2, 
  main = "A. Raw data",
  xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}

#Filter with filterByExpr function
keep <- filterByExpr(y, min.count = 10) #The filterByExpr function will only keep genes that have at least min.count counts per million (CPM) in at least n samples, where n is the number of samples in the smallest group (the pheromone group, in this case, where n = 6).
y_Filtered <- y[keep, keep.lib.sizes=F] #The keep.lib.sizes=F argument specifies that we want to recalculate the sample library sizes after removing the low-count genes from the data.

#Plot new logCPM reads after filtering
lcpm <- cpm(y_Filtered, log = T)
#Select this all loop and run it. if you run it line by line it won't work
plot(density(lcpm[, 1]),
     col = col[1],
     ylim = c(0, 0.6), 
     las = 2,
     main = "B. filterByExpr",
     xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples) {
  den <- density(lcpm[, i])
  lines(den$x, den$y, col = col[i])
}

#The next step is to normalize the data by RNA composition, in order to address problems caused by the counts in a sample being dominated by only a few highly expressed genes. This can cause other genes in the sample to exhibit artificially low count values. The normalization factors calculated using the calcNormFactors function are then used to scale the counts to avoid such problems.
#Calculate the normalization factors for each sample:
y_TMM <- calcNormFactors(y_Filtered, method="TMM")
#y_TMM <- calcNormFactors(y_Filtered, method="TMM")
#y_TMM <- calcNormFactors(y_Filtered, method="RLE")
#y_TMM <- calcNormFactors(y_Filtered, method="upperquartile")


#To view the normalization factors calculated for each sample (as well as library sizes), we can take a look at the samples attribute of our DGEList, y:
y_TMM$samples

#Let's check the effect
lcpm_y_Filtered <- cpm(y_Filtered, log = TRUE) #calculate lcpm of non-normalized data
lcpm_y_TMM <- cpm(y_TMM, log = TRUE) #calculate lcpm normalized data

#boxplot to check effect normalization (without outliers)
par(mfrow=c(2,2))
boxplot(lcpm_y_Filtered, las= 3, outline = FALSE, main = "Not normalized - No outliers", ylab = "Log-CPM")
boxplot(lcpm_y_TMM, las = 3, outline = FALSE, main = "Normalized - No outliers", ylab = "Log-CPM")

#boxplot to check effect normalization (with outliers)
boxplot(lcpm_y_Filtered, las= 3, outline = TRUE, main = "Not normalized - Outliers", ylab = "Log-CPM")
boxplot(lcpm_y_TMM, las = 3, outline = TRUE, main = "Normalized - Outliers", ylab = "Log-CPM")

#Dispersion estimation
#We can analyze dispersion estimation using pairwise contrasts, or using a generalized linear model. Here is demonstrated how to do pairwise.
#Estimate the dispersion across all genes:
#Using the grouping, we can now do an additional correction for the mean-variance correlation
y_MV <- estimateDisp(y_TMM)

y_MV$samples$lib.size

dev.off()
#We can also view the dispersion estimates using the BCV plot:
#The results of this plot give us an idea about the variance across all genes that are lowly expressed all the way to highly expressed. Normally, the more lowly expressed genes will have larger variation compared to the more highly expressed genes. (“Tagwise” is gene by gene dispersion, “common” is dispersion aross the whole dataset and “trend” takes groups of similarily expressed genes plot dispersion.)
plotBCV(y_MV)
dev.off()
