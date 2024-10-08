################################################
##   Source AskoR and set working directory   ##
##  (don't do if you use AskoR on Galaxy!)    ##
################################################
# Removes all objects from the current workspace (R memory)
rm(list=ls())
# source AskoR.R file
source("/directory/where/you/downloaded/the/file/AskoR.R")
# defined your workspace
setwd("/path/to/workspace/")

##############################################
##                Parameters                ##
##############################################
parameters<-Asko_start()

# Data and input files descriptions (if you are using AskoR in your R environment)
#--------------------------------------------------------------------------
# WARNING: All the input files must be in the same folder
#          called "input" (case sensitive)!
#--------------------------------------------------------------------------
parameters$analysis_name = "DEG_testPack"             # output directory name (default AskoRanalysis, do not put space!)
parameters$fileofcount = "CountsMatrix.txt"           # matrix of count for all samples/conditions
parameters$sep = "\t"                                 # field separator for count files or count matrix
parameters$annotation = "Genes_annotations.txt"       # file containing the functional annotations of each gene
parameters$geneID2GO_file = "GO_annotations.txt"      # GO annotation files
parameters$contrast_file = "Contrasts.txt"            # matrix of different contrasts desired
parameters$sample_file = "Samples_CountsMatrix.txt"   # file describing the samples
parameters$rm_sample = c("AC3R2","BC3R3")             # bad sample(s) !

# Specific Data import for use of AskoR on Galaxy
#--------------------------------------------------------------------------
# WARNING: All the input files must be in your Galaxy history!
#--------------------------------------------------------------------------
gx_get(1)                                             # Import dataset "1" from your Galaxy history
gx_get(2)                                             # Import dataset "2" from your Galaxy history
gx_get(3)                                             # Import dataset "3" from your Galaxy history
gx_get(4)                                             # Import dataset "4" from your Galaxy history
gx_get(5)                                             # Import dataset "5" from your Galaxy history

parameters$dir_path="./"                              # Workspace IF YOU ARE USINF AskoR ON GALAXY

parameters$analysis_name = "DEG_testPack"             # output directory name (default AskoRanalysis, do not put space!)
parameters$fileofcount = "1"                          # "1" if your matrix of count for all samples/conditions is your dataset "1" (if not, put the number of the dataset in your history)
parameters$sep = "\t"                                 # field separator for count files or count matrix
parameters$annotation = "2"                           # "2" if your file containing the functional annotations of each gene is your dataset "2" (if not, put the number of the dataset in your history)
parameters$geneID2GO_file = "3"                       # "3" if your GO annotation files is your dataset "3" (if not, put the number of the dataset in your history)
parameters$contrast_file = "4"                        # "4" if your matrix of different contrasts desired is your dataset "4" (if not, put the number of the dataset in your history)
parameters$sample_file = "5"                          # "5" if your file describing the samples is your dataset "5" (if not, put the number of the dataset in your history)
parameters$rm_sample = c("AC3R2","BC3R3")             # bad sample(s) !

# Options for data processing and their analyzes
#--------------------------------------------------------------------------
parameters$threshold_cpm = 0.5                        # CPM's threshold (default 0.5)
parameters$replicate_cpm = 3                          # Minimum number of replicates (default 3)
parameters$threshold_FDR = 0.05                       # FDR threshold (default 0.05)
parameters$threshold_logFC = 0                        # logFC threshold (default 0)  - if parameters$threshold_logFC != 0, AskoR will use glmTreat because use of "treat" is more accurate (see edgeR user guide section 2.13).
parameters$normal_method = "TMM"                      # normalization method (TMM/RLE/upperquartile/none) (default TMN)
parameters$p_adj_method = "fdr"                        # p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) (default fdr)
parameters$glm = "qlf"                                # GLM method (lrt/qlf) (default qlf)
parameters$CompleteHeatmap = FALSE                    # Generate Complete heatmap on all normalized genes (default FALSE, select "TRUE" only if your have less 30000 genes and enough RAM on your computer)
parameters$logFC = TRUE                               # logFC in the summary table (default TRUE)
parameters$FC = TRUE                                  # FC in the summary table (default TRUE)
parameters$logCPM = FALSE                             # logCPm in the summary table (default FALSE)
parameters$FDR = TRUE                                 # FDR in the summary table (default TRUE)
parameters$LR = FALSE                                  # LR in the summary table (default FALSE)
parameters$Sign = TRUE                                # Significance (1/0/-1) in the summary table (default TRUE)
parameters$Expression = TRUE                          # Significance expression in the summary table (default TRUE)
parameters$mean_counts = TRUE                         # Mean counts in the summary table (default TRUE)
parameters$norm_counts = FALSE                        # Generate files with mormalized counts (default FALSE)

# for legend of density plot
#-----------------------------------
parameters$densbotmar = 20                            # Set bottom margin of density plot to help position the legend (default 20)
parameters$densinset = 0.45                           # Set position the legend in bottom density graphe (default 0.45)
parameters$legendcol = 6                              # Set numbers of column for legends (default 6)

# Visualization of results from differential expression analyses
#-----------------------------------
parameters$plotMD = FALSE                              # Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
parameters$plotVO = FALSE                              # Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)
parameters$glimMD = FALSE                              # Glimma - Interactif Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
parameters$glimVO = FALSE                              # Glimma - Interactif Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)

########################################
##  Loading the data from the samples ##
########################################
##### load data #####
cat("\nRun DE analysis\n")
data<-loadData(parameters)
cat("\n\nChecking Data content:\n")
data$samples
data$contrast
data$design
head(data$dge$counts,n=4)

cat("Total number of genes : ", dim(data$dge$counts)[1], "\n")
cat("Total number of samples : ", dim(data$dge$counts)[2], "\n\n")
cat("summary of CPM by samples\n")
summary(edgeR::cpm(data$dge))
cat("\n")

##### asko files #####
asko_data<-asko3c(data, parameters)
cat("\nChecking Asko Data : condition, contrast, context.\n")
asko_data$condition ; cat("\n")
asko_data$contrast  ; cat("\n")
asko_data$context   ; cat("\n")

##### filtering #####
cat("\nFiltering genes with more than ", parameters$threshold_cpm, " CPM in ",parameters$replicate_cpm,"samples\n")
asko_filt<-GEfilt(data, parameters)
cat("Total number of filtered genes : ", dim(asko_filt$counts)[1], "\n\n")

##### normalization #####
asko_norm<-GEnorm(asko_filt, asko_data, data, parameters)

##### correlation #####
GEcorr(asko_norm,parameters)

##### DGE analysis #####
cat("\n\nDifferential expressions analysis\n")
resDEG<-DEanalysis(asko_norm, data, asko_data, parameters)

##### Venn diagram #####
# My list
parameters$compaVD = c("AC1vsAC2-AC1vsAC3-AC2vsAC3",
                       "BC1vsBC2-BC1vsBC3-BC2vsBC3",
                       "AC1vsBC1-AC2vsBC2-AC3vsBC3")

# graph type "all"
parameters$VD = "all"
VD(resDEG, parameters, asko_data)

# graph type "up"
parameters$VD = "up"
VD(resDEG, parameters, asko_data)

# graph type "down"
parameters$VD = "down"
VD(resDEG, parameters, asko_data)

# graph type "both"
parameters$compaVD = c("AC1vsBC1-AC2vsBC2",
                       "AC1vsBC1-AC3vsBC3",
                       "AC2vsBC2-AC3vsBC3")
parameters$VD = "both"
VD(resDEG, parameters, asko_data)

###### UpsetR Graphs #####
# My list
parameters$upset_list = c("AC1vsAC2-AC1vsAC3-AC2vsAC3",
                          "BC1vsBC2-BC1vsBC3-BC2vsBC3",
                          "AC1vsBC1-AC2vsBC2-AC3vsBC3")

# graphs type "all"
parameters$upset_basic = "all"
parameters$upset_type = "all"
UpSetGraph(resDEG, data, parameters)

# graphs type "mixed"
parameters$upset_basic = "mixed"
parameters$upset_type = "mixed"
UpSetGraph(resDEG, data, parameters)

# graphs type "up"
parameters$upset_basic = "up"
parameters$upset_type = "up"
UpSetGraph(resDEG, data, parameters)

# graphs type "down"
parameters$upset_basic = "down"
parameters$upset_type = "down"
UpSetGraph(resDEG, data, parameters)

##### Enrichment Analysis #####
# Parameters for enrichment
#----------------------------------------------------------------------
parameters$GO_threshold = 0.05               # the significant threshold used to filter p-values (default 0.05)
parameters$GO_min_num_genes = 10             # the minimum number of genes for each GO terms (default 10)
parameters$GO = "NULL"                       # gene set chosen for analysis 'up', 'down', 'both', or NULL (default NULL)
parameters$GO_algo = "weight01"              # algorithms for runTest function ("classic", "elim", "weight", "weight01", "lea", "parentchild") (default weight01)
parameters$GO_stats = "fisher"               # statistical tests for runTest function ("fisher", "ks", "t", "globaltest", "sum", "ks.ties") (default fisher)

# Parameters for visualization
#----------------------------------------------------------------------
parameters$Ratio_threshold = 0               # the min ratio for display GO in graph (default 0)
parameters$GO_max_top_terms = 10             # the maximum number of GO terms plot (default 10)
parameters$GO_min_sig_genes = 1              # the minimum number of significant gene(s) behind the enriched GO-term in graph (default 1)

# Run analysis on all contrasts
#----------------------------------------------------------------------
GOenrichment(resDEG, data, parameters)

# Run analysis on a gene list defined by the user
#----------------------------------------------------------------------
list=rownames(resDEG[1:1000,])              # contains a list of genes
GOenrichment(resDEG, data, parameters, list, "TitleOfTheList")


##### Co-Expression Analysis #####
# Parameters for gene clustering
#----------------------------------------------------------------------
parameters$coseq_data = "ExpressionProfiles"     # Perform clustering on transformed profiles based on normalized cpm counts (choose "LogScaledData" if you prefer to clusterize log2cpm counts and don't forget to set coseq_transformation to "none" in this case)
parameters$coseq_model = "kmeans"                # (default kmeans)
parameters$coseq_transformation = "clr"          # (default clr)
parameters$coseq_ClustersNb = 2:25               # (default : auto (select the best number automatically between 2 to 25))
parameters$coseq_HeatmapOrderSample = FALSE      # Choose TRUE if you prefer keeping your sample order than clusterizing samples in heatmap  (default FALSE)

# Parameters for for GO enrichment in clusters
#----------------------------------------------------------------------
parameters$GO_threshold = 0.05
parameters$GO_min_num_genes = 10
parameters$GO_algo = "weight01"
parameters$GO_stats = "fisher"

### Parameters for visualization
parameters$Ratio_threshold = 0
parameters$GO_max_top_terms = 10
parameters$GO_min_sig_genes = 0

# Run analysis on all DE genes (DE genes in 1 contrast at least)
#----------------------------------------------------------------------
clust<-ClustAndGO(asko_norm, resDEG, parameters, data)


# Include NON DE genes in an "artificial" cluster in the graphs produced by ClustAndGO
#----------------------------------------------------------------------
IncludeNonDEgenes_InClustering(data, asko_norm, resDEG, parameters, clust)
# Run analysis on a gene list defined by the user
#----------------------------------------------------------------------
list=rownames(resDEG[1:1000,])
ClustAndGO(asko_norm,resDEG,parameters, data, list, "TitleOfTheList")


##### Extract specific information on a gene list defined by the user #####
# -------------------------------------------------------------------------
list=rownames(resDEG[1:50,])

# Without clustering information
#----------------------------------------------------------------------
GeneInfo_OnList(list, resDEG, data, parameters, "TitleOfTheList")

# With clustering information (clust object produced by ClustAndGO function)
#----------------------------------------------------------------------
GeneInfo_OnList(list, resDEG, data, parameters, "TitleOfTheList",  clust)

# With selection of conditions and/or contrasts to be represented in the heatmap (clust object is not mandatory)
#----------------------------------------------------------------------
conditionsToDraw = c("AC1", "AC2", "AC3")                 # select conditions
contrastToDraw = c("AC1vsAC2","AC1vsAC3","AC2vsAC3")      # select contrasts
GeneInfo_OnList(list, resDEG, data, parameters, "TitleOfTheList", clust, contrasts=contrastToDraw)  # graph with selected contrasts and all conditions
GeneInfo_OnList(list, resDEG, data, parameters, "TitleOfTheList", clust, conditions=conditionsToDraw)  # graph with selected conditions and all contrasts
GeneInfo_OnList(list, resDEG, data, parameters, "TitleOfTheList", clust, conditions=conditionsToDraw, contrasts=contrastToDraw) ## graph with selected contrasts and selected conditions

