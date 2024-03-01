###################################################################################################################################################
###################################################################################################################################################
#######################################################                                     #######################################################
#######################################################  G E N E   S E T   A N A L Y S I S  ####################################################### 
#######################################################                                     #######################################################
###################################################################################################################################################
###################################################################################################################################################

#Extracting biological information from a volcano plot or from a result list with p-values and fold changes is quite challenging. 
#Even comparisons between strains that are quite similar can yield hundreds of differentially expressed genes. Manually examining all the genes would
#be quite time consuming, and the interpretation would be highly dependent on the scientist’s area of expertise. 
#Furthermore, gene-wise analysis fail to give a global overview of which cell components, processes, functions, or pathways are differentially 
#regulated between the cell-types being studied.
#To overcome these challenges, Gene Set Analysis (GSA) was developed. In a GSA, the gene sets are defined based on prior biological knowledge, 
#such as published data about biochemical pathways, biological processes, cellular function, and cellular components. 
#The goal of a GSA is to determines whether a defined gene set shows statistically significant differences between two condition in analysis. 
#Several methods for performing gene set analysis have been developed, and each of them relies on different assumptions. 
#In this workflow, I decided to perform GSA using the R package PIANO, developed by our division.

#Step 1. Make a list with the comparison you wanna make (in our case WTs vs the rest, including each other)
contrast.matrix.GSA <- makeContrasts(
  WT31vsCEN.PK = WT31 - CEN.PK,
  WT109vsCEN.PK = WT109 - CEN.PK,
  WT31vsKE6 = WT31 - KE6,
  WT109vsKE6 = WT109 - KE6,
  WT31vsEtOH = WT31 - EtOH,
  WT109vsEtOH = WT109 - EtOH,
  WT109vsWT31 = WT109 - WT31,
  levels = colnames(design.matrix))

vfit.GSA <- contrasts.fit(vfit, 
                          contrasts = contrast.matrix.GSA) 

vfit.GSA <- eBayes(vfit.GSA)

comparison_list.GSA <- list()
for (i in 1:ncol(contrast.matrix.GSA)) {
  name = colnames(contrast.matrix.GSA)[i]
  tmp = topTreat(vfit.GSA, coef = i, n = Inf) 
  comparison_list.GSA[[name]] <- tmp
}

#Step 2. Load the required functions - These were made by Simone and Ed to perform several kinds of GSAs in order to find the best/common results 
#among them. Loading this functions is basically the same as loading a new library.
source('DEanalysisMauri/functions/consGSA.R')
source('DEanalysisMauri/functions/consGSAplot.R')

#Step 3. Import gene set collection
#To perform a GSA, it is necessary to have a gene set collection (GSC). In our case it's the GSC of S. cerevisiae.
#In a GCS we can find four columns:
#Gene Sets: contains the id that uniquely identify a specific gene set.
#Genes: contains the gene names belonging to a specific gene set.
#Level: contains the levels to which a specific GO term is associated. It can have different levels such as biological process, molecular function, or cellular component.
#Term: contains the name of the GO term, which usually is quite descriptive.
#For well studied organisms such as S. cerevisiae, public databases with the gene sets exists (Saccharomyces genome database). 
#For less studied organisms it might be necessary to generate a gene set collection.
gene_sets <- read.csv2("DEanalysisMauri/Gene_Set.csv") #import the gene set collection

#After loading the file containing the GSC, we need to run the loadGSC function from the PIANO package to create a GSC object that will be 
#called by the consGSA function that runs the gene set analysis.
#We have three levels of analysis, and we want to perform a GSA for each
#We need to split the GSA by level
gene_sets <- list(geneSetBP = loadGSC(gene_sets %>%
                                        filter(Level %in% "biological_process") %>%
                                        select(c("genes", "Term")),
                                      type = "data.frame"),
                  geneSetMF = loadGSC(gene_sets %>%
                                        filter(Level %in% "molecular_function") %>%
                                        select(c("genes", "Term")),
                                      type = "data.frame"),
                  geneSetCC = loadGSC(gene_sets %>%
                                        filter(Level %in% "cellular_component") %>%
                                        select(c("genes", "Term")),
                                      type = "data.frame"))

#Perform GSA
#Several methods have been developed for running a GSA. Since each of these methods has different assumptions, different methods can give a 
#different interpretation of the data. To overcome this problem, it is possible to run a consensus gene set analysis (consGSA function). 
#The consGSA function runs a gene set analysis for each of the following GSA methods: mean, median, sum, stouffer, tailStrength, fisher, maxmean. 
#Further details can be found here: https://pubmed.ncbi.nlm.nih.gov/23444143/
#The next snippet of code loops the consensus gene set analysis through all the comparisons we have stored in the comparison_list. 
#This will generate a nested list called GSA_results. This list contains for each comparison one list, which in turn contains the results of 
#the gene set analysis for biological process, molecular function, and cellular component.

GSA_results.WT31vsCEN <- list() #Create list to store results of the consGSA function
GSA_results.WT31vsKE <- list() #Create list to store results of the consGSA function
GSA_results.WT31vsEtOH <- list() #Create list to store results of the consGSA function
GSA_results.WT109vsCEN <- list() #Create list to store results of the consGSA function
GSA_results.WT109vsKE <- list() #Create list to store results of the consGSA function
GSA_results.WT109vsEtOH <- list() #Create list to store results of the consGSA function
GSA_results.WT109vsWT31 <- list() #Create list to store results of the consGSA function

TopTable <- as.data.frame(comparison_list.GSA$WT31vsCEN.PK, col.names = "", row.names = NULL) #create df for GSA
GSA_results.WT31vsCEN[[i]] <- list(
  BP = consGSA(TopTable, gene_sets$geneSetBP), #run consensus GSA for biological process
  MF = consGSA(TopTable, gene_sets$geneSetMF), #run consensus GSA for molecular function
  CC = consGSA(TopTable, gene_sets$geneSetCC)) #run consensus GSA for cellular component

TopTable <- as.data.frame(comparison_list.GSA$WT31vsKE6, col.names = "", row.names = NULL) #create df for GSA
GSA_results.WT31vsKE[[i]] <- list(
  BP = consGSA(TopTable, gene_sets$geneSetBP), #run consensus GSA for biological process
  MF = consGSA(TopTable, gene_sets$geneSetMF), #run consensus GSA for molecular function
  CC = consGSA(TopTable, gene_sets$geneSetCC)) #run consensus GSA for cellular component

TopTable <- as.data.frame(comparison_list.GSA$WT31vsEtOH, col.names = "", row.names = NULL) #create df for GSA
GSA_results.WT31vsEtOH[[i]] <- list(
  BP = consGSA(TopTable, gene_sets$geneSetBP), #run consensus GSA for biological process
  MF = consGSA(TopTable, gene_sets$geneSetMF), #run consensus GSA for molecular function
  CC = consGSA(TopTable, gene_sets$geneSetCC)) #run consensus GSA for cellular component

TopTable <- as.data.frame(comparison_list.GSA$WT109vsCEN.PK, col.names = "", row.names = NULL) #create df for GSA
GSA_results.WT109vsCEN[[i]] <- list(
  BP = consGSA(TopTable, gene_sets$geneSetBP), #run consensus GSA for biological process
  MF = consGSA(TopTable, gene_sets$geneSetMF), #run consensus GSA for molecular function
  CC = consGSA(TopTable, gene_sets$geneSetCC)) #run consensus GSA for cellular component

TopTable <- as.data.frame(comparison_list.GSA$WT109vsKE6, col.names = "", row.names = NULL) #create df for GSA
GSA_results.WT109vsKE[[i]] <- list(
  BP = consGSA(TopTable, gene_sets$geneSetBP), #run consensus GSA for biological process
  MF = consGSA(TopTable, gene_sets$geneSetMF), #run consensus GSA for molecular function
  CC = consGSA(TopTable, gene_sets$geneSetCC)) #run consensus GSA for cellular component

TopTable <- as.data.frame(comparison_list.GSA$WT109vsEtOH, col.names = "", row.names = NULL) #create df for GSA
GSA_results.WT109vsEtOH[[i]] <- list(
  BP = consGSA(TopTable, gene_sets$geneSetBP), #run consensus GSA for biological process
  MF = consGSA(TopTable, gene_sets$geneSetMF), #run consensus GSA for molecular function
  CC = consGSA(TopTable, gene_sets$geneSetCC)) #run consensus GSA for cellular component

TopTable <- as.data.frame(comparison_list.GSA$WT109vsWT31, col.names = "", row.names = NULL) #create df for GSA
GSA_results.WT109vsWT31[[i]] <- list(
  BP = consGSA(TopTable, gene_sets$geneSetBP), #run consensus GSA for biological process
  MF = consGSA(TopTable, gene_sets$geneSetMF), #run consensus GSA for molecular function
  CC = consGSA(TopTable, gene_sets$geneSetCC)) #run consensus GSA for cellular component

names(GSA_results.WT31vsCEN) <- c('WT31vsCEN.PK') #assign names to each element of the list
names(GSA_results.WT31vsKE) <- c('WT31vsKE6') #assign names to each element of the list
names(GSA_results.WT31vsEtOH) <- c('WT31vsEtOH') #assign names to each element of the list
names(GSA_results.WT109vsCEN) <- c('WT109vsCEN.PK') #assign names to each element of the list
names(GSA_results.WT109vsKE) <- c('WT109vsKE6') #assign names to each element of the list
names(GSA_results.WT109vsEtOH) <- c('WT109vsEtOH') #assign names to each element of the list
names(GSA_results.WT109vsWT31) <- c('WT109vsWT31') #assign names to each element of the list

############ plot #############

Pcutoff <- 0.01 #set p value cutoff
consGSA_plot_data.WT31vsCEN <- list() #create list in which to save the graphs
consGSA_plot_data.WT31vsKE <- list() #create list in which to save the graphs
consGSA_plot_data.WT31vsEtOH <- list() #create list in which to save the graphs
consGSA_plot_data.WT109vsCEN <- list() #create list in which to save the graphs
consGSA_plot_data.WT109vsKE <- list() #create list in which to save the graphs
consGSA_plot_data.WT109vsEtOH <- list() #create list in which to save the graphs
consGSA_plot_data.WT109vsWT31 <- list() #create list in which to save the graphs

#generate the graphs and save them in the list, along with the data to generate them
for (i in 1:length(GSA_results.WT31vsCEN)) { #loop through all the GSAs we have
  consGSA_plot_data.WT31vsCEN[[i]] <- list( #the results will be saved in a list
    BP = consGSAplot(resList = GSA_results.WT31vsCEN[[i]][[1]], #create the consensous plot for BP
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Biological Process - rankScore 20", " ,pval < ", Pcutoff)),
    
    MF = consGSAplot(resList =  GSA_results.WT31vsCEN[[i]][[2]], #create the consensous plot for MF
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Molecular Function - rankScore 20" , " ,pval < ", Pcutoff)),
    
    CC = consGSAplot(resList = GSA_results.WT31vsCEN[[i]][[3]], #create the consensous plot for CC
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Cellular Component - rankScore 20", " ,pval < ", Pcutoff))
  )
  
  names(consGSA_plot_data.WT31vsCEN)[i] <- names(GSA_results.WT31vsCEN)[i] #assign the names to the list of graphs
}

for (i in 1:length(GSA_results.WT31vsKE)) { #loop through all the GSAs we have
  consGSA_plot_data.WT31vsKE[[i]] <- list( #the results will be saved in a list
    BP = consGSAplot(resList = GSA_results.WT31vsKE[[i]][[1]], #create the consensous plot for BP
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Biological Process - rankScore 20", " ,pval < ", Pcutoff)),
    
    MF = consGSAplot(resList =  GSA_results.WT31vsKE[[i]][[2]], #create the consensous plot for MF
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Molecular Function - rankScore 20" , " ,pval < ", Pcutoff)),
    
    CC = consGSAplot(resList = GSA_results.WT31vsKE[[i]][[3]], #create the consensous plot for CC
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Cellular Component - rankScore 20", " ,pval < ", Pcutoff))
  )
  
  names(consGSA_plot_data.WT31vsKE)[i] <- names(GSA_results.WT31vsKE)[i] #assign the names to the list of graphs
}

for (i in 1:length(GSA_results.WT31vsEtOH)) { #loop through all the GSAs we have
  consGSA_plot_data.WT31vsEtOH[[i]] <- list( #the results will be saved in a list
    BP = consGSAplot(resList = GSA_results.WT31vsEtOH[[i]][[1]], #create the consensous plot for BP
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Biological Process - rankScore 20", " ,pval < ", Pcutoff)),
    
    MF = consGSAplot(resList =  GSA_results.WT31vsEtOH[[i]][[2]], #create the consensous plot for MF
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Molecular Function - rankScore 20" , " ,pval < ", Pcutoff)),
    
    CC = consGSAplot(resList = GSA_results.WT31vsEtOH[[i]][[3]], #create the consensous plot for CC
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Cellular Component - rankScore 20", " ,pval < ", Pcutoff))
  )
  
  
  names(consGSA_plot_data.WT31vsEtOH)[i] <- names(GSA_results.WT31vsEtOH)[i] #assign the names to the list of graphs
}

for (i in 1:length(GSA_results.WT109vsCEN)) { #loop through all the GSAs we have
  consGSA_plot_data.WT109vsCEN[[i]] <- list( #the results will be saved in a list
    BP = consGSAplot(resList = GSA_results.WT109vsCEN[[i]][[1]], #create the consensous plot for BP
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Biological Process - rankScore 20", " ,pval < ", Pcutoff)),
    
    MF = consGSAplot(resList =  GSA_results.WT109vsCEN[[i]][[2]], #create the consensous plot for MF
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Molecular Function - rankScore 20" , " ,pval < ", Pcutoff)),
    
    CC = consGSAplot(resList = GSA_results.WT109vsCEN[[i]][[3]], #create the consensous plot for CC
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Cellular Component - rankScore 20", " ,pval < ", Pcutoff))
  )
  
  names(consGSA_plot_data.WT109vsCEN)[i] <- names(GSA_results.WT109vsCEN)[i] #assign the names to the list of graphs
}

for (i in 1:length(GSA_results.WT109vsKE)) { #loop through all the GSAs we have
  consGSA_plot_data.WT109vsKE[[i]] <- list( #the results will be saved in a list
    BP = consGSAplot(resList = GSA_results.WT109vsKE[[i]][[1]], #create the consensous plot for BP
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Biological Process - rankScore 20", " ,pval < ", Pcutoff)),
    
    MF = consGSAplot(resList =  GSA_results.WT109vsKE[[i]][[2]], #create the consensous plot for MF
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Molecular Function - rankScore 20" , " ,pval < ", Pcutoff)),
    
    CC = consGSAplot(resList = GSA_results.WT109vsKE[[i]][[3]], #create the consensous plot for CC
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Cellular Component - rankScore 20", " ,pval < ", Pcutoff))
  )
  
  names(consGSA_plot_data.WT109vsKE)[i] <- names(GSA_results.WT109vsKE)[i] #assign the names to the list of graphs
}

for (i in 1:length(GSA_results.WT109vsEtOH)) { #loop through all the GSAs we have
  consGSA_plot_data.WT109vsEtOH[[i]] <- list( #the results will be saved in a list
    BP = consGSAplot(resList = GSA_results.WT109vsEtOH[[i]][[1]], #create the consensous plot for BP
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Biological Process - rankScore 20", " ,pval < ", Pcutoff)),
    
    MF = consGSAplot(resList =  GSA_results.WT109vsEtOH[[i]][[2]], #create the consensous plot for MF
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Molecular Function - rankScore 20" , " ,pval < ", Pcutoff)),
    
    CC = consGSAplot(resList = GSA_results.WT109vsEtOH[[i]][[3]], #create the consensous plot for CC
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Cellular Component - rankScore 20", " ,pval < ", Pcutoff))
  )
  
  names(consGSA_plot_data.WT109vsEtOH)[i] <- names(GSA_results.WT109vsEtOH)[i] #assign the names to the list of graphs
}

for (i in 1:length(GSA_results.WT109vsWT31)) { #loop through all the GSAs we have
  consGSA_plot_data.WT109vsWT31[[i]] <- list( #the results will be saved in a list
    BP = consGSAplot(resList = GSA_results.WT109vsWT31[[i]][[1]], #create the consensous plot for BP
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Biological Process - rankScore 20", " ,pval < ", Pcutoff)),
    
    MF = consGSAplot(resList =  GSA_results.WT109vsWT31[[i]][[2]], #create the consensous plot for MF
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Molecular Function - rankScore 20" , " ,pval < ", Pcutoff)),
    
    CC = consGSAplot(resList = GSA_results.WT109vsWT31[[i]][[3]], #create the consensous plot for CC
                     rankScore = 20, 
                     Pcutoff = Pcutoff, 
                     distinct = 'distinct', 
                     title = paste0("Cellular Component - rankScore 20", " ,pval < ", Pcutoff))
  )
  
  names(consGSA_plot_data.WT109vsWT31)[i] <- names(GSA_results.WT109vsWT31)[i] #assign the names to the list of graphs
}


#For each comparison arrange the three (BP, MF,CC) graphs and print the final plot
print(ggarrange(consGSA_plot_data.WT31vsCEN[[i]][["BP"]][["plot"]], 
                consGSA_plot_data.WT31vsCEN[[i]][["MF"]][["plot"]], 
                consGSA_plot_data.WT31vsCEN[[i]][["CC"]][["plot"]], 
                ncol=1, 
                nrow=3, 
                common.legend = TRUE, 
                legend = "bottom") %>%
        annotate_figure(top = text_grob('WT31 vs CEN.PK', 
                                        color = "red", face = "bold", size = 18)))


print(ggarrange(consGSA_plot_data.WT31vsKE[[i]][["BP"]][["plot"]], 
                consGSA_plot_data.WT31vsKE[[i]][["MF"]][["plot"]], 
                consGSA_plot_data.WT31vsKE[[i]][["CC"]][["plot"]], 
                ncol=1, 
                nrow=3, 
                common.legend = TRUE, 
                legend = "bottom") %>%
        annotate_figure(top = text_grob('WT31 vs KE6-12', 
                                        color = "red", face = "bold", size = 18)))


print(ggarrange(consGSA_plot_data.WT31vsEtOH[[i]][["BP"]][["plot"]], 
                consGSA_plot_data.WT31vsEtOH[[i]][["MF"]][["plot"]], 
                consGSA_plot_data.WT31vsEtOH[[i]][["CC"]][["plot"]], 
                ncol=1, 
                nrow=3, 
                common.legend = TRUE, 
                legend = "bottom") %>%
        annotate_figure(top = text_grob('WT31 vs EtOH red', 
                                        color = "red", face = "bold", size = 18)))


print(ggarrange(consGSA_plot_data.WT109vsCEN[[i]][["BP"]][["plot"]], 
                consGSA_plot_data.WT109vsCEN[[i]][["MF"]][["plot"]], 
                consGSA_plot_data.WT109vsCEN[[i]][["CC"]][["plot"]], 
                ncol=1, 
                nrow=3, 
                common.legend = TRUE, 
                legend = "bottom") %>%
        annotate_figure(top = text_grob('WT109 vs CEN.PK', 
                                        color = "red", face = "bold", size = 18)))


print(ggarrange(consGSA_plot_data.WT109vsKE[[i]][["BP"]][["plot"]], 
                consGSA_plot_data.WT109vsKE[[i]][["MF"]][["plot"]], 
                consGSA_plot_data.WT109vsKE[[i]][["CC"]][["plot"]], 
                ncol=1, 
                nrow=3, 
                common.legend = TRUE, 
                legend = "bottom") %>%
        annotate_figure(top = text_grob('WT109 vs KE6-12', 
                                        color = "red", face = "bold", size = 18)))


print(ggarrange(consGSA_plot_data.WT109vsEtOH[[i]][["BP"]][["plot"]], 
                consGSA_plot_data.WT109vsEtOH[[i]][["MF"]][["plot"]], 
                consGSA_plot_data.WT109vsEtOH[[i]][["CC"]][["plot"]], 
                ncol=1, 
                nrow=3, 
                common.legend = TRUE, 
                legend = "bottom") %>%
        annotate_figure(top = text_grob('WT109 vs EtOH red', 
                                        color = "red", face = "bold", size = 18)))


print(ggarrange(consGSA_plot_data.WT109vsWT31[[i]][["BP"]][["plot"]], 
                consGSA_plot_data.WT109vsWT31[[i]][["MF"]][["plot"]], 
                consGSA_plot_data.WT109vsWT31[[i]][["CC"]][["plot"]], 
                ncol=1, 
                nrow=3, 
                common.legend = TRUE, 
                legend = "bottom") %>%
        annotate_figure(top = text_grob('WT109 vs WT31', 
                                        color = "red", face = "bold", size = 18)))


#Export each dataframe contained in each list
write.csv(consGSA_plot_data.WT31vsCEN$WT31vsCEN.PK$BP$df,"DEanalysisMauri/WT31vsCEN.BP.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT31vsKE$WT31vsKE$BP$df,"DEanalysisMauri/WT31vsKE.BP.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT31vsEtOH$WT31vsEtOH$BP$df,"DEanalysisMauri/WT31vsEtOH.BP.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsCEN$WT109vsCEN.PK$BP$df,"DEanalysisMauri/WT109vsCEN.BP.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsKE$WT109vsKE$BP$df,"DEanalysisMauri/WT109vsKE.BP.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsEtOH$WT109vsEtOH$BP$df,"DEanalysisMauri/WT109vsEtOH.BP.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsWT31$WT109vsWT31$BP$df,"DEanalysisMauri/WT109vsWT31.BP.csv", row.names = FALSE)

write.csv(consGSA_plot_data.WT31vsCEN$WT31vsCEN.PK$MF$df,"DEanalysisMauri/WT31vsCEN.MF.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT31vsKE$WT31vsKE$MF$df,"DEanalysisMauri/WT31vsKE.MF.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT31vsEtOH$WT31vsEtOH$MF$df,"DEanalysisMauri/WT31vsEtOH.MF.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsCEN$WT109vsCEN.PK$MF$df,"DEanalysisMauri/WT109vsCEN.MF.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsKE$WT109vsKE$MF$df,"DEanalysisMauri/WT109vsKE.MF.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsEtOH$WT109vsEtOH$MF$df,"DEanalysisMauri/WT109vsEtOH.MF.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsWT31$WT109vsWT31$MF$df,"DEanalysisMauri/WT109vsWT31.MF.csv", row.names = FALSE)

write.csv(consGSA_plot_data.WT31vsCEN$WT31vsCEN.PK$CC$df,"DEanalysisMauri/WT31vsCEN.CC.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT31vsKE$WT31vsKE$CC$df,"DEanalysisMauri/WT31vsKE.CC.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT31vsEtOH$WT31vsEtOH$CC$df,"DEanalysisMauri/WT31vsEtOH.CC.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsCEN$WT109vsCEN.PK$CC$df,"DEanalysisMauri/WT109vsCEN.CC.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsKE$WT109vsKE$CC$df,"DEanalysisMauri/WT109vsKE.CC.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsEtOH$WT109vsEtOH$CC$df,"DEanalysisMauri/WT109vsEtOH.CC.csv", row.names = FALSE)
write.csv(consGSA_plot_data.WT109vsWT31$WT109vsWT31$CC$df,"DEanalysisMauri/WT109vsWT31.CC.csv", row.names = FALSE)
