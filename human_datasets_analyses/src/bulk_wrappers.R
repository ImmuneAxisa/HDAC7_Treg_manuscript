## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------

##########################
## Packages and options ##
##########################


library(knitr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
# library(RColorBrewer)
library(viridis)
library(scales)
library(ggrepel)
library(tximport)
library(edgeR)
library(limma)
library(DESeq2)
library(RUVSeq)
library(EnsDb.Hsapiens.v86)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(biomaRt)
library(Hmisc)
library(DEFormats)
library(statmod)
#library(Cairo)
library(vsn)
#library(kableExtra)
library(gplots)
library(cowplot)
library(ggExtra)
library(sva)
library(FactoMineR)
library(factoextra)
library(ggsci)




##################
## ggplot theme ##
##################

# From Matthew R Lincoln

theme_mrl <- function(x = 1) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(size = 12*x,face = "bold", angle = 45, vjust = 0.6),
      axis.text.y = element_text(size = 12*x,face = "bold"),
      axis.title.x = element_text(size = 12*x,face = "bold"),
      axis.title.y = element_text(size = 12*x,face = "bold"),
      strip.background = element_rect(fill="gray20", colour="gray20", linetype="solid"),
                                      strip.text = element_text(size=14*x, colour="white", face="bold"),
                                      legend.title = element_text(size=14*x, face = "bold"),
                                      legend.text = element_text(size=12*x, color="gray20", face="bold"),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
      plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
      plot.caption = element_text(hjust = 0, face = "italic")
                                      )
}


## ----gene type filtering functions, cache = F---------------------------------------------------------------------------------------------


##################################################
## Download transcript-gene models from Ensembl ##
##################################################
edb <- EnsDb.Hsapiens.v86
tx2gene <- genes(edb, return.type = "data.frame")

###########################################################
## Remove unwanted genes entries fron a DGEList object ####
###########################################################

# Remove genes on the X chromosome. Works with DGEList and dds objects
RemoveX <- function(dgeList, tx2gene = tx2gene) {
  genes_to_remove <- tx2gene %>%
    dplyr::filter(seq_name == 'X') %>%
    dplyr::select("gene_id") %>%
    unname() %>%
    unlist()
  
  return(dgeList[!rownames(dgeList) %in% genes_to_remove, ])
}

# Remove genes on the Y chromosome. Works with DGEList and dds objects
RemoveY <- function(dgeList, tx2gene = tx2gene) {
  genes_to_remove <- genes(edb, return.type = "data.frame") %>%
    dplyr::filter(seq_name == 'Y') %>%
    dplyr::select("gene_id") %>%
    unname() %>%
    unlist()
  
  return(dgeList[!rownames(dgeList) %in% genes_to_remove, ])
}


# Keep only protein_coding and lincRNA biotypes. It will remove a series of non-coding RNAs, uncharaterized transcripts and TCR&IG segments.
# Works with DGEList and dds objects
OnlyCodingAndlinc <- function(dgeList, tx2gene = tx2gene) {
  genes_to_remove <- tx2gene %>%
    dplyr::filter(gene_biotype != "protein_coding" & gene_biotype != "lincRNA") %>%
    dplyr::select("gene_id") %>%
    unname() %>%
    unlist()
  
  return(dgeList[!rownames(dgeList) %in% genes_to_remove, ])
}


#Remove TCR and Ig segments
RemoveIgTcr <- function(dgeList, tx2gene = tx2gene) {
  genes_to_remove <- tx2gene %>%
    dplyr::filter(grepl("TR_", gene_biotype) | grepl("IG_", gene_biotype)) %>%
    dplyr::select("gene_id") %>%
    unname() %>%
    unlist()
  
  return(dgeList[!rownames(dgeList) %in% genes_to_remove, ])
}




## ----DGEList clean function, cache = F----------------------------------------------------------------------------------------------------


cleanDGEList <- function(DGEList, MIN_INDIVIDUALS = 3, MIN_CPM = 10, normMethod = "TMM", ggPsize = 1){
  
  ############################################Filterting low expression reads ##########################################
  # Log-CPM transform unfiltereted data:
  DGEList.prefiltered <- as.data.frame(cpm(DGEList, log = TRUE))

  DGEList.prefiltered <- DGEList.prefiltered %>%
    rownames_to_column(var = "Gene") %>%
    mutate(Stage = "Pre-filtering")

  # Filter reads
  DGEList.filtered <- DGEList[rowSums(cpm(DGEList) > MIN_CPM) >= MIN_INDIVIDUALS,, keep.lib.sizes=FALSE]
  DGE.filtered <- as.data.frame(cpm(DGEList.filtered, log = TRUE)) %>%
    rownames_to_column(var = "Gene") %>%
    mutate(Stage = "Post-filtering")
  
  # Build matrix of expression densities for plotting:
  expression.density <- bind_rows(DGEList.prefiltered, DGE.filtered) %>%
    gather(Sample, logCPM, -Gene, -Stage) %>%
    mutate(Stage = factor(Stage, levels = c("Pre-filtering", "Post-filtering")))
  
  #plot before/after filtering
  p1 <- ggplot(expression.density, aes(logCPM)) +
          geom_density(aes(group = Sample, color = Sample)) +
          geom_vline(xintercept = 0, linetype = 2, colour = "gray50") +
          scale_colour_viridis(discrete = TRUE) +
          labs(x = "log(CPM)", y = "Density") +
          facet_grid(. ~ Stage) +
          theme_mrl(ggPsize) +
          guides(colour = guide_legend(nrow = 4)) +
          theme(legend.position = "none") +
          ggtitle(paste0("Filtering of gene counts: Genes with low counts are removed. \n Include genes that are present at >", MIN_CPM," counts per million \n in at least ", MIN_INDIVIDUALS, " samples"))
  
  ##################################### Library normalization ##########################################
  
  
  # Prepare log-CPM
lcpm.pre <- as.data.frame(cpm(DGEList.filtered, log = TRUE)) %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Stage = "Pre-normalization")

# Normalize:
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = normMethod)

lcpm.post <- as.data.frame(cpm(DGEList.filtered.norm, log = TRUE)) %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Stage = "Post-normalization")


norm.data <- bind_rows(lcpm.pre, lcpm.post) %>%
  gather(Sample, logCPM, -Gene, -Stage) %>%
  mutate(Stage = factor(Stage, levels = c("Pre-normalization", "Post-normalization")))


#Plot before/after normalization
p2 <- ggplot(norm.data, aes(x = Stage, color = Stage)) +
              geom_boxplot(aes(y = logCPM), fill = "grey20", outlier.colour = "grey50", outlier.size = 0.5) +
              facet_wrap(~ Sample, ncol = 18) +
              labs(y = "log(CPM)") +
              theme_mrl(ggPsize*0.6) +
              theme(axis.text.x = element_blank()) +
              theme(legend.position = "bottom") +
              theme(panel.spacing=unit(0.001, "lines")) +
              ggtitle(paste0("Normalization of libraries by ", normMethod))

p3 <- ggplot(norm.data, aes(x = Sample, color = Sample)) +
              geom_boxplot(aes(y = logCPM), fill = "grey20", outlier.colour = "grey50", outlier.size = 0.5) +
              facet_grid(. ~ Stage) +
              labs(y = "log(CPM)") +
              theme_mrl(ggPsize) +
              theme(axis.text.x = element_blank()) +
              theme(legend.position = "none") +
              ggtitle(paste0("Normalization of libraries by ", normMethod))
  
  print(p1)
  print(p2)
  print(p3)
  return(DGEList.filtered.norm)
}



## ----Limma functions, cache = F-----------------------------------------------------------------------------------------------------------

###########################################
## obtain differentially expressed genes ##
###########################################

deg.limma.voom <- function(input.dgelist, design, coef = NULL, contrast = NULL, lfc = 0, 
                           TREAT = FALSE, trend = FALSE, robust = FALSE, tx2gene = NULL, top =50) {
  # To estimate confounders, perform initial association without covariates and obtain least
  # significant genes:
  # design <- model.matrix(~ phenotype, data=input.dgelist$samples)
  voom <- voom(input.dgelist, design, plot = TRUE, save.plot = TRUE)
  # voom <- voom(input.dgelist$counts, design, plot = TRUE, save.plot = TRUE)
  fit <- lmFit(voom, design)
  if (TREAT){
    fit <- treat(fit, lfc = lfc, trend = trend, robust = robust)
    deg <- topTreat(fit, coef = coef, sort.by = "P", n = Inf)
  } else {
    fit <- eBayes(fit, trend = trend, robust = robust)
    deg <- topTable(fit, coef = coef, sort.by = "P", n = Inf)
  }
  plotSA(fit, main = "Final model: mean-variance trend")
  plotMD(fit, coef = coef, status = deg$adj.P.val*100)
  # Obtain differentially expressed genes:  
  if (is.null(tx2gene) == FALSE) {
    deg <- as.data.frame(deg) %>% 
      tibble::rownames_to_column("gene_id") %>% 
      left_join(tx2gene %>% dplyr::select(gene_id, gene_name, gene_biotype), by = "gene_id")
  }
  print(kable(top_n(deg, n = -top, wt = P.Value) , 
              caption = "Limma top Differentially expressed genes", 
              longtable = TRUE, 
              booktabs = TRUE)) 
  return(deg)
}

##############################################################################
## obtain differentially expressed genes after modelling unwanted variation ##
##############################################################################

deg.limma.voom.ruv <- function(input.dgelist, design, ruv.confounder.panel.size = 1000, top = 50) {
  # To estimate confounders, perform initial association without covariates and obtain least
  # significant genes:
  # design <- model.matrix(~ phenotype, data=input.dgelist$samples)
  voom <- voom(input.dgelist, design, plot = FALSE)
  fit <- lmFit(voom, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  least.sig.genes <- topTable(fit, coef = "conditionhdac7WT", sort.by = "P", n = Inf) %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    arrange(-P.Value) %>%
    dplyr::slice(1:ruv.confounder.panel.size) %>%
    dplyr::select(gene_id) %>%
    unlist() %>% unname()
  
  # Voom transform and calculate residuals:
  # voom.ruv <- voom(input.dgelist$counts, design, plot = FALSE, save.plot = TRUE)
  # res <- sapply(1:nrow(voom.ruv$E), function(x) lm(voom.ruv$E[x,] ~ input.dgelist$samples$phenotype, weights = voom.ruv$weights[x,])$residuals)
  
  # Use residuals and least significant genes to obtain W matrix from RUV:
  set <- RUVg(voom$E, rownames(voom$E) %in% least.sig.genes, k = 1, isLog = TRUE)
  # set <- RUVr(voom.ruv$E, rownames(voom.ruv$E) %in% least.sig.genes, k=1, t(res), isLog=TRUE)
  W <- set$W
  
  # Use W matrix to model unwanted variation:
  par(mfrow = c(1, 2))
  input.dgelist$samples$W = W
  # design.model = model.matrix(~ phenotype + W, data=input.dgelist$samples)
  design.model <- as.data.frame(design)
  design.model$W <- W
  
  voom.model <- voom(input.dgelist$counts, design.model, plot = TRUE, save.plot = TRUE)
  fit.model <- lmFit(voom.model, design)
  fit.model <- eBayes(fit.model, trend = TRUE, robust = TRUE)
  plotSA(fit.model, main = "Final model: mean-variance trend")

  # Obtain differentially expressed genes:  
  deg <- topTable(fit.model, coef = "conditionhdac7WT", sort.by = "P", n = Inf) %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    left_join(read_tsv("data/Homo_sapiens.GRCh38.88.primary_assembly.gene_symbols.txt"), by = "gene_id") %>%
    dplyr::select(gene_name, gene_id, everything())
  
  print(kable(top_n(deg, n = -top, wt = adj.P.Val) , caption = "top Differentially expressed genes", longtable = TRUE, booktabs = TRUE))
  
  return(deg)
}





## ----EdgeR function, cache = F------------------------------------------------------------------------------------------------------------
DEGedgeR <- function(DGEList, designMat, coef, contrast = NULL, dispersion = "common", 
                     lfc = 0, TREAT = FALSE, tx2gene = NULL, top = 50) {
  #Variance plot for each sample
  plotMDS.DGEList(DGEList, 
                  main = "Multidimensional scaling plot of distances between digital gene expression profiles",
                  cex.main = 0.8)
  #Estimate dispesion with common or trended model
  if (dispersion == "common") {
    DGEList <- estimateGLMCommonDisp(DGEList, design=designMat)
  }
  if (dispersion == "trended") {
    DGEList <- estimateGLMTrendedDisp(DGEList, design=designMat)
  } 
  DGEList <- estimateGLMTagwiseDisp(DGEList, design=designMat)
  plotBCV(DGEList, main = "BCV estimate")
  #fitting curve
  fit <- glmFit(DGEList, designMat)
  p <- gof(fit, plot =  TRUE)
  
  #DGE testing LRT vs TREAT
  if (TREAT) {
    results <- glmTreat(fit, coef = coef, lfc = lfc)
  }
  else {
    results <- glmLRT(fit, coef = coef)
    #results$table <- results$table %>% 
     # tibble::rownames_to_column("gene") %>% 
      #filter(abs(logFC)>lfc) %>% 
      #tibble::column_to_rownames("gene")
  }
  #Calculate FDR and output DEG table
  SortedResults <- topTags(results, adjust.method = "fdr", n = nrow(results$table))
  plotSmear(results, de.tags =  rownames(topTags(results, adjust.method = "fdr", n = 1000, p.value = 0.1)$table), cex = 1)
  #If a tx2gene correspondance table for the gene_id is provided, join it to the results table
  if (is.null(tx2gene) == FALSE) {
    degs <- SortedResults$table %>% 
      tibble::rownames_to_column("gene_id") %>% 
      left_join(tx2gene %>% dplyr::select(gene_id, gene_name, gene_biotype), by = "gene_id") 
  }
  else {
    degs <- as.data.frame(SortedResults$table)
  }
  print(qqplot)
  print(kable(top_n(degs, n = -top, wt = PValue), 
              caption = "EdgeR top Differentially expressed genes", 
              longtable = TRUE, booktabs = TRUE)) 
  return(degs)
}


## ----DEGS with DESeq with dds object input, cache = F-------------------------------------------------------------------------------------


degDESeq <- function(dds, Min_count = 1, alpha = 0.1, contrast = NULL, fitType = "local", tx2gene = tx2gene, top = 50, lfcT = 0) {
  #remove low counts genes
  dds <- dds[ rowSums(counts(dds)) >= Min_count, ]
  #Get DEGs
  dds <- DESeq(dds, fitType = fitType)
  #plot dispersion estimates
  plotDispEsts(dds)
  cat("\n\n")
  #Gen-wise standard deviation after fitting
  meanSdPlot(assay(vst(dds, blind = FALSE)))
  cat("\n\n")
  #Extract DEGs based on contrast and lfc threshold
  degs <- results(dds, alpha=alpha, contrast=contrast, lfcThreshold = lfcT)
  #Calculate shrunk LFC with the "normal" default DESeq2 method.
  shrink <-lfcShrink(dds, alpha = alpha, contrast=contrast, type = "normal", lfcThreshold = lfcT)
  #Plot MA plots with summary metrics.
  metrics <- capture.output(summary(degs))
  DESeq2::plotMA(degs, ylim = c(-5.1,5.1), 
                 main = paste0("Raw LFC:  ",metrics[4],"    ;    ", metrics[5]), 
                 sub = paste0(metrics[6],"    ;    ", metrics[7]), 
                 cex =0.4,
                 cex.main =0.8,
                 cex.sub =0.8)
  cat("\n\n")
  DESeq2::plotMA(shrink, ylim = c(-2.1,2.1), 
                 main = paste0("Shrunken LFC:  ",metrics[4],"    ;    ", metrics[5]), 
                 sub = paste0(metrics[6],"    ;    ", metrics[7]), 
                 cex =0.4,
                 cex.main =0.8,
                 cex.sub =0.8)
  cat("\n\n")
  #format shrink object for merging
  shrink <- as.data.frame(shrink) %>%
      tibble::rownames_to_column("gene_id") %>% 
      mutate(Slfc = log2FoldChange) %>%
      dplyr::select(gene_id, Slfc)
  #merge raw and shrunk data
  degs <- as.data.frame(degs) %>%
      tibble::rownames_to_column("gene_id") %>% 
      left_join(tx2gene %>% dplyr::select(gene_id, gene_name, gene_biotype), by = "gene_id") %>%
      left_join(shrink, by = "gene_id") %>%
      arrange(padj)
  #kable(top_n(degs, n = -top, wt = padj) , 
              #caption = "DESeq2 top Differentially expressed genes", 
              #longtable = TRUE, 
              #booktabs = TRUE)
  #Normalized expression distribution
  p<- ggplot() +
    geom_point(alpha = 0.5, data = degs %>% dplyr::filter(padj >=0.1), aes(x = log10(baseMean), y = lfcSE, color = padj < 0.1)) +
    geom_point(alpha = 0.5, data = degs %>% dplyr::filter(padj <0.1), aes(x = log10(baseMean), y = lfcSE, color = padj < 0.1)) +
    theme(legend.position = "bottom") +
    ggsci::scale_color_jco()
  print(plot_grid(ggMarginal(p, fill = "grey")))
  cat("\n\n")
  #PLot lfcSE vs LFC
  p<- ggplot() +
    geom_point(alpha = 0.1, data = degs %>% dplyr::filter(padj >=0.1), aes(x = lfcSE, y = abs(log2FoldChange), color = padj < 0.5)) +
    geom_point(alpha = 0.1, data = degs %>% dplyr::filter(padj <0.1), aes(x = lfcSE, y = abs(log2FoldChange), color = padj < 0.5)) +
    theme(legend.position = "bottom") +
    ggsci::scale_color_jco()
  print(plot_grid(ggMarginal(p, fill = "grey")))
  cat("\n\n")
return(degs)
}


degDESeqASHR <- function(dds, Min_count = 1, alpha = 0.1, contrast = NULL, fitType = "local", tx2gene = tx2gene, top = 50) {
  #remove low counts genes
  dds <- dds[ rowSums(counts(dds)) >= Min_count, ]
  #Get DEGs
  dds <- DESeq(dds, fitType = fitType)
  plotDispEsts(dds)
  
  meanSdPlot(assay(vst(dds, blind = FALSE)))
  
  degs <- results(dds, alpha=alpha, contrast=contrast)
  shrink <-lfcShrink(dds, alpha = alpha, contrast=contrast, type = "ashr", svalue = T)
  metrics <- capture.output(summary(degs))
  DESeq2::plotMA(degs, ylim = c(-5.1,5.1), 
                 main = paste0("Raw LFC:  ",metrics[4],"    ;    ", metrics[5]), 
                 sub = paste0(metrics[6],"    ;    ", metrics[7]), 
                 cex =0.4,
                 cex.main =0.8,
                 cex.sub =0.8)
  
  DESeq2::plotMA(shrink, ylim = c(-2.1,2.1), 
                 main = paste0("Shrunken LFC:  ",metrics[4],"    ;    ", metrics[5]), 
                 sub = paste0(metrics[6],"    ;    ", metrics[7]), 
                 cex =0.4,
                 cex.main =0.8,
                 cex.sub =0.8)
  
  shrink <- as.data.frame(shrink)
  colnames(shrink) <- paste0(colnames(shrink), "_ashr")
  shrink <- shrink %>%
      tibble::rownames_to_column("gene_id")
  degs <- as.data.frame(degs) %>%
      tibble::rownames_to_column("gene_id") %>% 
      left_join(tx2gene %>% dplyr::select(gene_id, gene_name, gene_biotype), by = "gene_id") %>%
      left_join(shrink, by = "gene_id") %>%
      arrange(padj)
  #kable(top_n(degs, n = -top, wt = padj) , 
              #caption = "DESeq2 top Differentially expressed genes", 
              #longtable = TRUE, 
              #booktabs = TRUE)
  #plot normalized expression distribution
  hist(log10(degs$baseMean), 
       xlab = "log10 base mean", main = "Base mean distribution", 
       breaks = 30)
  abline(v=1, col = "blue")
  #PLot lfcSE vs LFC
  p<- ggplot() +
    geom_point(alpha = 1, data = degDESeqVariant %>% dplyr::filter(padj >=0.5), aes(x = lfcSE, y = abs(log2FoldChange), color = padj < 0.5)) +
    geom_point(alpha = 0.5, data = degDESeqVariant %>% dplyr::filter(padj <0.5), aes(x = lfcSE, y = abs(log2FoldChange), color = padj < 0.5)) +
    theme(legend.position = "bottom") +
    ggsci::scale_color_jco()
  plot_grid(ggMarginal(p, fill = "grey"))
  
return(degs)
}





## ----volvano plot function, cache = F-----------------------------------------------------------------------------------------------------

plot.volcano <- function(results) {

### Selects results columns needed for the plot, log-transforms the p-values:
data <- data.frame(gene_id = results$gene_id,
                   gene_name = results$gene_name,
                   pvalue = -log10(results$pvalue), 
                   padj = results$padj,
                   lfc = results$Slfc)
# Removes NAs
data <- na.omit(data)

### Creates a seperate column FDR_cat, categorizing genes by FDR:
data <- data %>%
  mutate(FDR_cat = ifelse(data$padj < 0.05, "FDR < 0.05", ifelse(data$padj < 0.1, "FDR < 0.1","FDR > 0.1"))) 

### Graphical parameters (ggplot2):

vplot <- ggplot(data, aes(x = lfc, y = pvalue)) + # sets general plot structure
  geom_point(aes(color = factor(FDR_cat)), size = 1.7, alpha = 0.5, na.rm = T) + # tweaks the look of the dots
  theme_minimal(base_size = 16) + # sets visual theme of the plot (background, fonts etc.)
  theme(legend.title=element_text(size = 0),legend.text=element_text(size = 10)) + # legend formatting
  xlab(expression(log[2]("Fold Change"))) + # x-axis label
  ylab(expression(-log[10]("p-value"))) + # y-axis label
  geom_vline(xintercept = 0, colour = "black") + # sets an x axis intercept at 0
  geom_hline(yintercept = 1.3, colour = "black") + # adds a y axis intercept at 1.3 (p = 0.05) 
  scale_color_manual(values = c("FDR < 0.05" = "red", 
                                "FDR < 0.1" = "orange", 
                                "FDR > 0.1" = "dark grey")) # changes dot color based on FDR group
  


### Selects 10 genes with lowest FDR:
top_labelled <- top_n(data, n = 10, wt = (padj*-1))

### Labels dots corresponiding to 10 lowest FDR genes with gene names:
vplot<-vplot + geom_text_repel(data = top_labelled, 
                          mapping = aes(label = gene_name), 
                          size = 2,
                          fontface = 'bold', 
                          color = 'black',
                          box.padding = unit(0.5, "lines"),
                          point.padding = unit(0.5, "lines"))
    
return(vplot) # returns plot
}



## ----annotated plot function--------------------------------------------------------------------------------------------------------------

# function to plot DEGs and annotate selected genes with their gene names.
# takes as input the differential expression df output from DESeq2:
# Required columns: X, Y, FDR_cat, gene_name
# a character vector or genes to annotate,
# the x and y limits of the plots.
annotated.plot <- function(dataset, annotatedGenes = c(), S = 1, bins = 20, nudgeX=0, nudgeY=0) {

diffplot <- ggplot(dataset, aes(x = X, y = Y)) + # sets general plot structure
  geom_point(data = subset(dataset, FDR_cat %in% c("FDR > 0.1")), 
             size = 1.7, alpha = 0.5, na.rm = T, shape = 21, fill = "dark grey") +
  stat_density2d(aes(color = stat(level)), bins = bins) +
  geom_point(data = subset(dataset, FDR_cat %in% c("FDR < 0.1")) , 
             size = 2, alpha = 0.5, na.rm = T, shape = 21, fill = "orange") +
  theme(legend.title=element_text(size = 10),legend.text=element_text(size = 10), legend.position = "none") # legend formatting
  #geom_smooth( method = "lm", colour = "red") +
  #theme(aspect.ratio = 1)
  #geom_abline(slope = 1, intercept = 0.5, linetype = "dotted") +
  #geom_abline(slope = 1, intercept = -0.5, linetype = "dotted") +

diffplot1 <- diffplot + 
  theme_mrl(S) + # sets visual theme of the plot (background, fonts etc.)
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6), legend.position = "none")
### Selects 10 genes with lowest FDR:
#top_labelled <- top_n(CombinedDEGsplots, n = 50, wt = (padjIT*-1))
top_labelled <- dataset %>%
  dplyr::filter(gene_name %in% annotatedGenes)

### Labels dots corresponiding to 50 lowest FDR genes with gene names:
diffplot2<-diffplot + geom_label_repel(data = top_labelled, 
                          mapping = aes(label = gene_name), 
                          size = 2.5,
                          #force = 1,
                          #ylim = c(-1,1.5,1.5),
                          #xlim = c(-1.5,1.5,1),
                          nudge_x= nudgeX*top_labelled$X,
                          nudge_y= nudgeY*top_labelled$Y,
                          max.iter = 5000,
                          fontface = 'bold', 
                          color = 'black',
                          segment.colour = "black",
                          box.padding = unit(0.95, "lines"),
                          label.padding = unit(0.15, "lines"),
                          point.padding = unit(0.25, "lines")) +
  ggnewscale::new_scale_fill() +
  geom_point(data = top_labelled, aes(fill = factor(FDR_cat)), shape = 21, color = "red", alpha = 1, stroke = 1) +
  scale_fill_manual(values = c("FDR < 0.1" = "orange", 
                                "FDR > 0.1" = "dark grey"))+# changes dot color based on FDR group
                    
cat("\n\n")

diffplot3 <- diffplot2 + 
  theme_mrl(S) + # sets visual theme of the plot (background, fonts etc.)
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6), legend.position = "none")

return(list(diffplot1, diffplot3, plot_grid(diffplot1, diffplot3, ncol =2)))
}







#### utils


#function to turn named list into a data frame where 
#one column contains the name of the vector
#and the other is the vector values, one by row
#from https://gist.github.com/aammd/9ae2f5cce9afd799bafb
list_to_df <- function(listfordf){
  if(!is.list(listfordf)) stop("it should be a list")
  
  df <- list(list.element = listfordf)
  class(df) <- c("tbl_df", "data.frame")
  attr(df, "row.names") <- .set_row_names(length(listfordf))
  
  if (!is.null(names(listfordf))) {
    df$name <- names(listfordf)
  }
  
  df
}


#this function uses biomartR to convert mouse ensembl ids to human and get the human gene name.
#the input needs to be an atomic vector. Ouput is a dataframe, with mouse id, human id and human gene name.
convertMouseGeneList <- function(x) {
  
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://jul2018.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "http://jul2018.archive.ensembl.org")
  
  
  genesV2 = getLDS(
    attributes = c("ensembl_gene_id"),
    filters = "ensembl_gene_id",
    values = x,
    mart = mouse,
    attributesL = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    martL = human,
    uniqueRows = T
  )
  #humanx <- unique(genesV2[, 2])
  humanx <- genesV2
  
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  return(humanx)
}


