#####################################################################################################################################
#####################################################################################################################################
#################################################                                   #################################################
#################################################      D G E    A N A L Y S I S     ################################################# 
#################################################                                   #################################################
#####################################################################################################################################
#####################################################################################################################################


#CREATE DESIGN AND CONTRAST MATRICES
#To find relationships that are of biological interest, appropriate statistical models need to be applied to the data. Differential expression analysis 
#through limma uses linear models to determine the size and direction of the changes in gene expression. The modelling process through limma requires the
#use of a design matrix that defines the experiment structure and variables. After sub-setting the data, I use only the “Strain” variable to create a
#design matrix, since I am only interested in the effect of this factor. I also decided to use a means model (~ 0 + group) instead of a mean-reference 
#model. As documented by Law, C.W. et al. (2020) the two models are equivalent for categorical variables (e.g. cell-type).
#Create model matrix
strain <- y_MV$samples$group #y_MV equals to "Subset" df in Simone's original script.
design.matrix <- model.matrix(~ 0 + strain, data = y_MV) #define design matrix
rownames(design.matrix) <- rownames(y_MV$samples)
colnames(design.matrix) <- gsub("strain", "", colnames(design.matrix))

#Now we need to create a contrast matrix, that will be later used to specify which comparisons between group levels are to be extracted from the model.
#To build a constrast matrix we can use the makeContrasts function. Note that the second factor in the contrast matrix, CEN.PK in our first case, will be the reference. 
#This means that a gene will be upregulated/downregulated in the first factor compared to CEN.PK.
#Create contrast matrix
contrast.matrix.vsCEN <- makeContrasts(
  KE6vsCEN.PK = KE6 - CEN.PK,
  EtOHvsCEN.PK = EtOH - CEN.PK,
  WT31vsCEN.PK = WT31 - CEN.PK,
  WT109vsCEN.PK = WT109 - CEN.PK,
  levels = colnames(design.matrix))

contrast.matrix.vsKE <- makeContrasts(
  EtOHvsKE6 = EtOH - KE6,
  WT31vsKE6 = WT31 - KE6,
  WT109vsKE6 = WT109 - KE6,
  levels = colnames(design.matrix))

contrast.matrix.vsEtOH <- makeContrasts(
  WT31vsEtOH = WT31 - EtOH,
  WT109vsEtOH = WT109 - EtOH,
  levels = colnames(design.matrix))

contrast.matrix.vsWT31 <- makeContrasts(
  WT109vsWT31 = WT109 - WT31,
  levels = colnames(design.matrix))

#VOOM TRANSFORMATION
#Limma works with voom-transformed data, so before fitting a model to the data we need to perform a voom transformation.
#A voom transformation consists of the following steps:

#1. Counts are transformed to log2-CPM, where the “per million reads” is defined based on the normalization factors that were calculated with the calcNormFactors function.
#2. A linear model is fitted to the log2-CPM for each gene (gene-wise linear model fitting). The model borrows information from the design matrix and takes into account experimental design, experimental conditions, and replicates.
#3. When fitting the model, a standard deviation of the residuals for each gene is calculated. The residual is the difference between the observed value and the estimated value (for example, the mean).
#4. A robust trend is then fitted to the residual standard deviations as a function of the average log2-CPM for each gene (red line in the voom-plot). The standard deviation trend is then used to calculate weights for each observation. Typically, the “voom-plot” shows a decreasing trend between the means and variances. Each observation with its log2-CPM values and associated weights is the input into limma’s differential expression pipeline.
#perform voom transformation
v <- voom(y_MV, 
          design.matrix,
          plot = TRUE) 

#FIT THE MODEL
#Linear modelling in limma is carried out using the lmFit and contrasts.fit functions. These functions were originally written for application to 
#microarrays, but after performing a voom transformation they can be used for RNA-seq data.
#The lmFit function fits a linear model using weighted least squares for each gene. Weighted least squares is a generalization of ordinary least 
#squares in which the variance of observations, that was estimated by the voom transformation, is incorporated into the regression.
#The contrasts.fit function allows us to extract the fold-change between each gene in the groups specified by the contrast.matrix.
vfit <- lmFit(v, 
              design.matrix)

vfit.vsCEN <- contrasts.fit(vfit, 
                            contrasts = contrast.matrix.vsCEN) 

vfit.vsKE <- contrasts.fit(vfit, 
                           contrasts = contrast.matrix.vsKE) 

vfit.vsEtOH <- contrasts.fit(vfit, 
                             contrasts = contrast.matrix.vsEtOH) 

vfit.vsWT31 <- contrasts.fit(vfit, 
                             contrasts = contrast.matrix.vsWT31)

#Next, empirical Bayes moderation is carried out by borrowing information across all genes to obtain more precise estimates of gene-wise variability. 
#The eBayes function “takes the linear model as input and computes moderated t-statistics, moderated F-statistic, and log-odds of differential 
#expression by empirical Bayes moderation of the standard errors towards a global value”.
vfit.vsCEN <- eBayes(vfit.vsCEN)
vfit.vsKE <- eBayes(vfit.vsKE)
vfit.vsEtOH <- eBayes(vfit.vsEtOH)
vfit.vsWT31 <- eBayes(vfit.vsWT31)


#For a quick look at differential expression levels, the number of significantly up- and down-regulated genes can be summarized in a table. 
#The default significance cutoff is based on adjusted p value and is set at 5%.
summary(decideTests(vfit.vsCEN, adjust.method = "BH", p.value = 0.01, lfc = 0)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsKE, adjust.method = "BH", p.value = 0.01, lfc = 0)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsEtOH, adjust.method = "BH", p.value = 0.01, lfc = 0)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsWT31, adjust.method = "BH", p.value = 0.01, lfc = 0)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.

summary(decideTests(vfit.vsCEN, adjust.method = "BH", p.value = 0.01, lfc = 1)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsKE, adjust.method = "BH", p.value = 0.01, lfc = 1)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsEtOH, adjust.method = "BH", p.value = 0.01, lfc = 1)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsWT31, adjust.method = "BH", p.value = 0.01, lfc = 1)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.

summary(decideTests(vfit.vsCEN, adjust.method = "BH", p.value = 0.01, lfc = 3.3)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsKE, adjust.method = "BH", p.value = 0.01, lfc = 3.3)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsEtOH, adjust.method = "BH", p.value = 0.01, lfc = 3.3)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.
summary(decideTests(vfit.vsWT31, adjust.method = "BH", p.value = 0.01, lfc = 3.3)) #Save each summary in a excel CSV file; the second obj added in the vfit.XXX element is the control, whereas the first obj is the one which is up-/down-regulated compared to the control. in this case we have DE6-12 as control and CEN.PK as sample. therefore the up-/down-regulation is of CEN.PK compared to KE and should be flipped when annotated.

#Display number of DE genes for each couple in a barplot
#First we do lfc=0, then lfc=1
myData <- read.csv(file.choose(), header = TRUE, sep = ";", dec = ",")#file named "DE genes counts in anaerobic conditions - Simone script" in the Mauri RNAseq analysis folder, containing the number of genes up-/down- regulated for each strain when compared to each control

myData$DE <- factor(myData$DE, levels = c("up","down"))
myData$Strain <- factor(myData$Strain, levels = c("KE6-12","Ethanol red",'LBCM31','LBCM109',"Ethanol red.",'LBCM31.','LBCM109.','LBCM31:','LBCM109:','LBCM109!'))

#lcf=0
plot1 <- ggplot(subset(myData, lfc == "0"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs CEN.PK', subtitle='p.value 0.01; lfc 0')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=30, hjust=1, vjust=1),
        axis.text.y = element_text(size=30),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot1


plot2 <- ggplot(subset(myData, vs == "KE" & lfc == "0"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs KE', subtitle='p.value 0.01; lfc 0')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot3 <- ggplot(subset(myData, vs == "Et" & lfc == "0"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs EtOH', subtitle='p.value 0.01; lfc 0')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot4 <- ggplot(subset(myData, vs == "WT31" & lfc == "0"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs WT31', subtitle='p.value 0.01; lfc 0')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15, face='bold'))

grid.arrange(plot1, plot2, plot3, plot4, ncol=4)

#lcf=1 (one)
plot1 <- ggplot(subset(myData, lfc == "1"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs CEN.PK', subtitle='p.value 0.01; lfc 0')+
  ylim(0,3500) +
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=30, hjust=1, vjust=1),
        axis.text.y = element_text(size=30),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot1

plot2 <- ggplot(subset(myData, vs == "KE" & lfc == "1"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs KE', subtitle='p.value 0.01; lfc 1')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot3 <- ggplot(subset(myData, vs == "Et" & lfc == "1"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs EtOH', subtitle='p.value 0.01; lfc 1')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot4 <- ggplot(subset(myData, vs == "WT31" & lfc == "1"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs WT31', subtitle='p.value 0.01; lfc 1')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15, face='bold'))

grid.arrange(plot1, plot2, plot3, plot4, ncol=4)

#lcf=3.3
plot1 <- ggplot(subset(myData, vs == "CENPK" & lfc == "3.3"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs CEN.PK', subtitle='p.value 0.01; lfc 3.3')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot2 <- ggplot(subset(myData, vs == "KE" & lfc == "3.3"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs KE', subtitle='p.value 0.01; lfc 3.3')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot3 <- ggplot(subset(myData, vs == "Et" & lfc == "3.3"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs EtOH', subtitle='p.value 0.01; lfc 3.3')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position='none',
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2))

plot4 <- ggplot(subset(myData, vs == "WT31" & lfc == "3.3"), aes(x=Strain, y=Mean, fill=DE)) + 
  geom_col() +
  geom_text(aes(label=Mean), position=position_stack(vjust=.5), size=5) +
  ggtitle('DE genes vs WT31', subtitle='p.value 0.01; lfc 3.3')+
  ylab("Number of genes") +
  xlab("Strain")+
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"), panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = 'black', size = 1.5), 
        axis.text.x = element_text(angle=45, size=10, hjust=1, vjust=1),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15, vjust=2),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15, face='bold'))

grid.arrange(plot1, plot2, plot3, plot4, ncol=4)

#MOVE RESULTS IN A LIST
#Since usually multiple comparisons are performed, it is convenient to move all the results in a list, in which each element is the result of one comparison. This will also make possible to easily loop downstream analysis to all the comparisons, without writing additional code. From now on, all the code will automatically loop though all the comparisons and perform the analysis on each of them automatically, without having to modify or write additional code.
#move the results of the comparisons in a list
comparison_list.vsCEN <- list()
for (i in 1:ncol(contrast.matrix.vsCEN)) {
  name = colnames(contrast.matrix.vsCEN)[i]
  tmp = topTreat(vfit.vsCEN, coef = i, n = Inf) 
  comparison_list.vsCEN[[name]] <- tmp
}

comparison_list.vsKE <- list()
for (i in 1:ncol(contrast.matrix.vsKE)) {
  name = colnames(contrast.matrix.vsKE)[i]
  tmp = topTreat(vfit.vsKE, coef = i, n = Inf) 
  comparison_list.vsKE[[name]] <- tmp
}

comparison_list.vsEtOH <- list()
for (i in 1:ncol(contrast.matrix.vsEtOH)) {
  name = colnames(contrast.matrix.vsEtOH)[i]
  tmp = topTreat(vfit.vsEtOH, coef = i, n = Inf) 
  comparison_list.vsEtOH[[name]] <- tmp
}

comparison_list.vsWT31 <- list()
for (i in 1:ncol(contrast.matrix.vsWT31)) {
  name = colnames(contrast.matrix.vsWT31)[i]
  tmp = topTreat(vfit.vsWT31, coef = i, n = Inf) 
  comparison_list.vsWT31[[name]] <- tmp
}

#In each list, for each comparison, we obtain a table. The columns of these tables contain:
#logFC: log2 fold change
#AveExpr: average expression of the gene in all the samples, expressed in log2-CPM
#t: logFC divided by its standard error
#p.value: raw p-value derived from testing that logFC is different from 0
#adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
#B: log-odds that gene is DE (arguably less useful than the other columns)

#Export each dataframe contained in each list
for (i in seq_along(comparison_list.vsCEN)) {
  filename = paste(i, ".csv")
  write.csv(comparison_list.vsCEN[[i]], filename)
}

for (i in seq_along(comparison_list.vsKE)) {
  filename = paste(i, ".csv")
  write.csv(comparison_list.vsKE[[i]], filename)
}

for (i in seq_along(comparison_list.vsEtOH)) {
  filename = paste(i, ".csv")
  write.csv(comparison_list.vsEtOH[[i]], filename)
}

for (i in seq_along(comparison_list.vsWT31)) {
  filename = paste(i, ".csv")
  write.csv(comparison_list.vsWT31[[i]], filename)
}

###################################################
######### E N H A N C E D   V O L C A N O #########
###################################################
#Volcano plots with p.value 0.01 and logFC 3.3
#Volcano plot vs CEN.PK 
remotes::install_github('kevinblighe/EnhancedVolcano', force=TRUE) #use it if you have troubles with fitting all the DEG in the plot


volc1.10updown <- EnhancedVolcano(comparison_list.vsCEN$KE6vsCEN.PK,
                                  lab = rownames(comparison_list.vsCEN$KE6vsCEN.PK),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('FDH1','GAL80','YPL276W',
                                                'AAD3','RDS1','HPF1','SPG1','GRE1','GAL2', 'GAL10', 'MAL12','GPP2','YAR028W', 
                                                'PRM7','REE1','HXT3','IMD2','STD1','ICS2','IMA1'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'KE vs CEN.PK',
                                  subtitle = "lfc referred to KE compared to CEN.PK",
                                  caption = "",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')

volc2.10updown <- EnhancedVolcano(comparison_list.vsCEN$EtOHvsCEN.PK,
                                  lab = rownames(comparison_list.vsCEN$EtOHvsCEN.PK),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('HXT5','YGR201C','SPG1',
                                                'MSC1','HPF1','FMP45','HBT1','ECM4','MBR1', 'YIP3', 'HOM3','TRP3','REE1', 
                                                'PRM7','YDL073W','YAR028W','ICY2','MAL11','ALD5','MAL12'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'EtOH vs CEN.PK',
                                  subtitle = "lfc referred to EtOH compared to CEN.PK",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')

volc3.10updown <- EnhancedVolcano(comparison_list.vsCEN$WT31vsCEN.PK,
                                  lab = rownames(comparison_list.vsCEN$WT31vsCEN.PK),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('HXT5','GRE1','SPG1',
                                                'ATG36','GTT1','HPF1','YGR201C','AYT1','FMP45', 'AAD3','MAL12','YAR028W','REE1','HOM3',
                                                'YER188W','PRM7','GPP2','MAL11','YIL014C-A','ICS2'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'WT31 vs CEN.PK',
                                  subtitle = "lfc referred to WT31 compared to CEN.PK",
                                  caption = "",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')

volc4.10updown <- EnhancedVolcano(comparison_list.vsCEN$WT109vsCEN.PK,
                                  lab = rownames(comparison_list.vsCEN$WT109vsCEN.PK),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('SPG1','YGR201C','MBR1',
                                                'MSC1','GTT1','GRE1','HPF1','SPS100','YDR018C', 'COS4','HOM3','IMA1','MAL12','ICY2','TRP3',
                                                'HXT3','YAR028W','GPP2','PRM7','STE2'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'WT31 vs CEN.PK',
                                  subtitle = "lfc referred to WT109 compared to CEN.PK",
                                  caption = "fold-of-change cutoff = 10; Adj.p-value cutoff = 0.01",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')


grid.arrange(volc1.10updown, volc2.10updown, volc3.10updown, volc4.10updown, ncol=2)

#Volcano plot vs KE6-12
volc5.10updown <- EnhancedVolcano(comparison_list.vsKE$EtOHvsKE,
                                  lab = rownames(comparison_list.vsKE$EtOHvsKE),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('YIP3','HXT5','PNS1',
                                                'COS1','DIA1','MPC3','MSC1','PIL1','LAP2', 'SSE2','SMF1','YPL276W','ZRT1','EXG1',
                                                'TRP3','KTR2','YPL067C','HIS3','AIF1','SSK22'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'EtOH vs KE',
                                  subtitle = "lfc referred to EtOH compared to KE",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')

volc6.10updown <- EnhancedVolcano(comparison_list.vsKE$WT31vsKE,
                                  lab = rownames(comparison_list.vsKE$WT31vsKE),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('HXT5','YFL054C','MAL13',
                                                'SPS100','YGR109W-B','YFL052W','PEX18','YHR210C','MAL11', 'FEX2','FDH1','GAL80','YPL276W','ZRT1',
                                                'YHB1','KTR2','RDS1','PHR1','YER186C','YNL058C'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'WT31 vs KE',
                                  subtitle = "lfc referred to WT31 compared to KE",
                                  caption = "",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')

volc7.10updown <- EnhancedVolcano(comparison_list.vsKE$WT109vsKE,
                                  lab = rownames(comparison_list.vsKE$WT109vsKE),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('SPS100','MSC1','YFL054C',
                                                'STF2',"YFL052W",'ACO1','SOL4','TKL2','POM33','YHR210C','FDH1','TMT1','YPL276W',
                                                'ZRT1','TRP3','HIS3','EXG1','ARG4','SSK22','PHR1'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'WT109 vs KE',
                                  subtitle = "lfc referred to WT109 compared to KE",
                                  caption = "fold-of-change cutoff = 10; Adj.p-value cutoff = 0.01",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')


grid.arrange(volc5.10updown, volc6.10updown, volc7.10updown, ncol=2)

#Volcano plot vs EtOH and WT31vsWT109
volc8.10updown <- EnhancedVolcano(comparison_list.vsEtOH$WT31vsEtOH,
                                  lab = rownames(comparison_list.vsEtOH$WT31vsEtOH),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('SMF1','IFM1','FLO10','MAL11',
                                                'YGL262W','SEA4','THI22','AIF1','PPR1','NFT1','YIP3','YER188W','DIA1','RRI2','SCW4','SSE2',
                                                'FDH1','ENA1','COS1','BAP2'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'WT31 vs EtOH',
                                  subtitle = "lfc referred to WT31 compared to EtOH",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')

volc9.10updown <- EnhancedVolcano(comparison_list.vsEtOH$WT109vsEtOH,
                                  lab = rownames(comparison_list.vsEtOH$WT109vsEtOH),
                                  x = 'logFC',
                                  y = 'adj.P.Val',
                                  selectLab = c('SMF1','QDR2','CTH1',
                                                'SUR2','AIF1','HSP33','YCR100C','YNR073C','DIA3','YGL262W','YIP3','BAP2','DIA1','FAT3','YER188W','HXT4','ENA1',
                                                'FDH1','COS1','IMA1'),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylim=c(0,25),
                                  xlim=c(-20,20),
                                  pCutoff = 0.01,
                                  FCcutoff = 3.3,
                                  title = 'WT109 vs EtOH',
                                  subtitle = "lfc referred to WT109 compared to EtOH",
                                  caption = "",
                                  pointSize = 4.0,
                                  maxoverlapsConnectors = Inf,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  legendPosition = 'none',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black')

volc10.10updown <- EnhancedVolcano(comparison_list.vsWT31$WT109vsWT31,
                                   lab = rownames(comparison_list.vsWT31$WT109vsWT31),
                                   x = 'logFC',
                                   y = 'adj.P.Val',
                                   selectLab = c('DOG2','BNA6','YHB1',
                                                 'SCW4','YPR015C','YER188W','YOL159C-A','RDS1','BNA5','SUR7','IMA1','MAL13','MAL11','UIP3','YLR460C',
                                                 'REE1','YAR029W','HXT7','FLO1','YIR042C'),
                                   xlab = bquote(~Log[2]~ 'fold change'),
                                   ylim=c(0,25),
                                   xlim=c(-20,20),
                                   pCutoff = 0.01,
                                   FCcutoff = 3.3,
                                   title = 'WT109 vs WT31',
                                   subtitle = "lfc referred to WT109 compared to WT31",
                                   caption = "fold-of-change cutoff = 10; Adj.p-value cutoff = 10e-6",
                                   pointSize = 4.0,
                                   maxoverlapsConnectors = Inf,
                                   labSize = 4.0,
                                   labCol = 'black',
                                   legendPosition = 'none',
                                   legendLabSize = 14,
                                   legendIconSize = 4.0,
                                   drawConnectors = TRUE,
                                   widthConnectors = 1.0,
                                   colConnectors = 'black')

grid.arrange(volc8.10updown, volc9.10updown, volc10.10updown, ncol=2)

#Volcano plots with p.value 0.01 and logFC 1
#Volcano plot vs CEN.PK 
volc1.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsCEN$KE6vsCEN.PK,
                                       lab = rownames(comparison_list.vsCEN$KE6vsCEN.PK),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('FDH1','GAL80','YPL276W',
                                                     'AAD3','RDS1','HPF1','SPG1','GRE1','GAL2', 'ECM4', 'MAL12','GPP2','YAR028W', 
                                                     'PRM7','REE1','HXT3','IMD2','STD1','ICS2','IMA1'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'KE vs CEN.PK',
                                       subtitle = "lfc referred to KE compared to CEN.PK",
                                       maxoverlapsConnectors = Inf,
                                       caption = "",
                                       pointSize = 4.0,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

volc2.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsCEN$EtOHvsCEN.PK,
                                       lab = rownames(comparison_list.vsCEN$EtOHvsCEN.PK),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('HXT5','YGR201C','SPG1',
                                                     'MSC1','HPF1','FMP45','HBT1','ECM4','MBR1', 'YIP3', 'HOM3','TRP3','REE1', 
                                                     'PRM7','YDL073W','YAR028W','DFG16','ICY2','MAL11','ALD5'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'EtOH vs CEN.PK',
                                       subtitle = "lfc referred to EtOH compared to CEN.PK",
                                       maxoverlapsConnectors = Inf,
                                       pointSize = 4.0,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

volc3.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsCEN$WT31vsCEN.PK,
                                       lab = rownames(comparison_list.vsCEN$WT31vsCEN.PK),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('HXT5','GRE1','SPG1',
                                                     'ATG36','GTT1','HPF1','YGR201C','ECM4', 'MBR1','AYT1','MAL12','YAR028W','REE1','HOM3',
                                                     'YER188W','PRM7','GPP2','MAL11','YIL014C-A','ICS2'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'WT31 vs CEN.PK',
                                       subtitle = "lfc referred to WT31 compared to CEN.PK",
                                       caption = "",
                                       maxoverlapsConnectors = Inf,
                                       pointSize = 4.0,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

volc4.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsCEN$WT109vsCEN.PK,
                                       lab = rownames(comparison_list.vsCEN$WT109vsCEN.PK),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('SPG1','YGR201C','MBR1',
                                                     'MSC1','GTT1','GRE1','HPF1','SPS100','YDR018C', 'COS4','HOM3','IMA1','MAL12','ICY2','TRP3',
                                                     'HXT3','YAR028W','GPP2','PRM7','STE2'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'WT31 vs CEN.PK',
                                       subtitle = "lfc referred to WT109 compared to CEN.PK",
                                       maxoverlapsConnectors = Inf,
                                       caption = "fold-of-change cutoff = 2; Adj.p-value cutoff = 0.01",
                                       pointSize = 4.0,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')


grid.arrange(volc1.10updown.lFC1, volc2.10updown.lFC1, volc3.10updown.lFC1, volc4.10updown.lFC1, ncol=2)

#Volcano plot vs KE6-12
volc5.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsKE$EtOHvsKE,
                                       lab = rownames(comparison_list.vsKE$EtOHvsKE),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('YIP3','HXT5','PNS1',
                                                     'COS1','DIA1','MPC3','MSC1','PIL1','LAP2', 'SSE2','SMF1','YPL276W','FKS1','ZRT1',
                                                     'EXG1','TRP3','KTR2','SLP1','YPL067C','MET22'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'EtOH vs KE',
                                       subtitle = "lfc referred to EtOH compared to KE",
                                       pointSize = 4.0,
                                       maxoverlapsConnectors = Inf,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

volc6.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsKE$WT31vsKE,
                                       lab = rownames(comparison_list.vsKE$WT31vsKE),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('HXT5','YFL054C','MAL13',
                                                     'YHR140W','SPS100','YGR109W-B','YFL052W','PEX18','YHR210C', 'MAL11','FDH1','GAL80','YPL276W','ZRT1',
                                                     'YHB1','AAD3','KTR2','RDS1','RIO1','PHR1'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'WT31 vs KE',
                                       subtitle = "lfc referred to WT31 compared to KE",
                                       caption = "",
                                       pointSize = 4.0,
                                       maxoverlapsConnectors = Inf,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

volc7.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsKE$WT109vsKE,
                                       lab = rownames(comparison_list.vsKE$WT109vsKE),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('SPS100','MSC1','YFL054C',
                                                     'STF2','YHR140W',"YFL052W",'ACO1','SOL4','TKL2','POM33','FDH1','TMT1','YPL276W',
                                                     'ZRT1','TRP3','HOM3','ARG7','HIS3','RDS1','EXG1'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'WT109 vs KE',
                                       subtitle = "lfc referred to WT109 compared to KE",
                                       caption = "fold-of-change cutoff = 2; Adj.p-value cutoff = 0.01",
                                       pointSize = 4.0,
                                       maxoverlapsConnectors = Inf,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')


grid.arrange(volc5.10updown.lFC1, volc6.10updown.lFC1, volc7.10updown.lFC1, ncol=2)

#Volcano plot vs EtOH and WT31vsWT109
volc8.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsEtOH$WT31vsEtOH,
                                       lab = rownames(comparison_list.vsEtOH$WT31vsEtOH),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('SMF1','GPI10','ATG36',
                                                     'RQC1','YDR061W','IFM1','YSP1','QDR2','PMP1','CYS4','YIP3','YER188W','DIA1','RRI2','SCW4','SSE2',
                                                     'FDH1','ENA1','COS1','BAP2'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'WT31 vs EtOH',
                                       subtitle = "lfc referred to WT31 compared to EtOH",
                                       pointSize = 4.0,
                                       maxoverlapsConnectors = Inf,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

volc9.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsEtOH$WT109vsEtOH,
                                       lab = rownames(comparison_list.vsEtOH$WT109vsEtOH),
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       selectLab = c('SMF1','QDR2','CTH1','ECL1','LST7','SUR2',
                                                     'COX10','PSR1','PMP1','PRM4','YIP3','BAP2','DIA1','FAT3','YER188W','MAL33','RRI2',
                                                     'HXT4','ENA1','FDH1'),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       ylim=c(0,25),
                                       xlim=c(-20,20),
                                       pCutoff = 0.01,
                                       FCcutoff = 1,
                                       title = 'WT109 vs EtOH',
                                       subtitle = "lfc referred to WT109 compared to EtOH",
                                       caption = "",
                                       pointSize = 4.0,
                                       maxoverlapsConnectors = Inf,
                                       labSize = 4.0,
                                       labCol = 'black',
                                       legendPosition = 'none',
                                       legendLabSize = 14,
                                       legendIconSize = 4.0,
                                       drawConnectors = TRUE,
                                       widthConnectors = 1.0,
                                       colConnectors = 'black')

volc10.10updown.lFC1 <- EnhancedVolcano(comparison_list.vsWT31$WT109vsWT31,
                                        lab = rownames(comparison_list.vsWT31$WT109vsWT31),
                                        x = 'logFC',
                                        y = 'adj.P.Val',
                                        selectLab = c('DOG2','GAL80','TOP1',
                                                      'BNA6','YHB1','COX10','UPS3','SCW4','YLR049C','SOP4','IMA1','MAL13','TMT1','STB4','MAL11',
                                                      'RPS25B','SER33','FAT3','SSM4','RIB3'),
                                        xlab = bquote(~Log[2]~ 'fold change'),
                                        ylim=c(0,25),
                                        xlim=c(-20,20),
                                        pCutoff = 0.01,
                                        FCcutoff = 1,
                                        title = 'WT109 vs WT31',
                                        subtitle = "lfc referred to WT109 compared to WT31",
                                        caption = "fold-of-change cutoff = 2; Adj.p-value cutoff = 0.01",
                                        pointSize = 4.0,
                                        maxoverlapsConnectors = Inf,
                                        labSize = 4.0,
                                        labCol = 'black',
                                        legendPosition = 'none',
                                        legendLabSize = 14,
                                        legendIconSize = 4.0,
                                        drawConnectors = TRUE,
                                        widthConnectors = 1.0,
                                        colConnectors = 'black')

grid.arrange(volc8.10updown.lFC1, volc9.10updown.lFC1, volc10.10updown.lFC1, ncol=2)

#find the DE genes up-/down-regulated in more than 1 strain
myData2 <- read.csv(file.choose(), header = TRUE, sep = ";", dec = ",") # csv. file "All top10 up-downreguöated genes - FC-2 or FC10 - Anaerobic" in "Mauri RNAseq analysis" folder

topDE.1 <- data.frame(table(myData2$up.CEN.10))
topDE.2 <- data.frame(table(myData2$down.CEN.10))
topDE.3 <- data.frame(table(myData2$up.KE.10))
topDE.4 <- data.frame(table(myData2$down.KE.10))
topDE.5 <- data.frame(table(myData2$up.EtOH.10))
topDE.6 <- data.frame(table(myData2$down.EtOH.10))
topDE.7 <- data.frame(table(myData2$up.WT31.10))
topDE.8 <- data.frame(table(myData2$down.WT31.10))
topDE.9 <- data.frame(table(myData2$up.WT109.10))
topDE.10 <- data.frame(table(myData2$down.WT109.10))

topDE.11 <- data.frame(table(myData2$up.CEN.2))
topDE.12 <- data.frame(table(myData2$down.CEN.2))
topDE.13 <- data.frame(table(myData2$up.KE.2))
topDE.14 <- data.frame(table(myData2$down.KE.2))
topDE.15 <- data.frame(table(myData2$up.EtOH.2))
topDE.16 <- data.frame(table(myData2$down.EtOH.2))
topDE.17 <- data.frame(table(myData2$up.WT31.2))
topDE.18 <- data.frame(table(myData2$down.WT31.2))
topDE.19 <- data.frame(table(myData2$up.WT109.2))
topDE.20 <- data.frame(table(myData2$down.WT109.2))

mylist <- list(up.CEN.10 = topDE.1, down.CEN.10 = topDE.2,
               up.KE.10 = topDE.3, down.KE.10 = topDE.4,
               up.EtOH.10 = topDE.5, down.EtOH.10 = topDE.6,
               up.WT31.10 = topDE.7, down.WT31.10 = topDE.8,
               up.WT109.10 = topDE.9, down.WT109.10 = topDE.10,
               up.CEN.2 = topDE.11, down.CEN.2 = topDE.12,
               up.KE.2 = topDE.13, down.KE.2 = topDE.14,
               up.EtOH.2 = topDE.15, down.EtOH.2 = topDE.16,
               up.WT31.2 = topDE.17, down.WT31.2 = topDE.18,
               up.WT109.2 = topDE.19, down.WT109.2 = topDE.20)

#Export each dataframe contained in each list
#it saves them in the documents named with numbers ("1", "2", etc...)
for (i in seq_along(mylist)) {
  filename = paste(i, ".csv")
  write.csv(mylist[[i]], filename)
}
