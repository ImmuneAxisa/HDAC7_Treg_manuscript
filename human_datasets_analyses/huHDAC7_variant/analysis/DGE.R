#####################################################################################################################
############################################## Set up ###############################################################
#####################################################################################################################

# Import packages and custom wrapper functions
source("src/bulk_wrappers.R")

# other libraries
library(rhdf5)



#####################################################################################################################
#################################### Import quantitation results ####################################################
#####################################################################################################################


# read in sample table
all.samples <- read_csv("huHDAC7_variant/data/GEO_table_STM_variant.csv") %>% spread(fastq_read, file)

#Convert all relevant factors for the linear regression to factors if not the case
all.samples$Batch <- as.factor(all.samples$Batch)
all.samples$Treatment <- factor(all.samples$Treatment, levels = c("ctl", "hdac7WT", "hdac7Mut"))
all.samples$sampleID <- as.factor(all.samples$sampleID)


sample.table.RNA <- all.samples %>%
  dplyr::filter(!is.na(RNA_FastQ1) & !is.na(RNA_FastQ2)) %>%
  mutate(RNA_Htseq_File = paste0("huHDAC7_variant/results/rna/quantitation/htseq-count/", Sample_Name, ".gtf"),
         RNA_RSEM_File = paste0("huHDAC7_variant/results/rna/quantitation/rsem/", Sample_Name, ".genes.results"))

# Import RSEM-quantitated reads with tximport:
txi.rsem <- tximport(sample.table.RNA$RNA_RSEM_File,
                     type = "rsem", countsFromAbundance = "lengthScaledTPM")
#prevents error when importing in DESeq
txi.rsem$length[txi.rsem$length == 0] <- 1

#Import data from rsem quantitation into a dds object
TR.dds <- DESeqDataSetFromTximport(txi.rsem, all.samples, design = ~1)

#Change the design after removing the deprecated levels of cell type and treatment factors.
#Mandatory to input the design only at the end otherwise the design matrix will be of the wrong size.
TR.dds@design <- ~Batch+Treatment


#Remove irrelevant gene entries: genes on the Y chromosome and Ig/TCR segments
TR.dds.curated <- RemoveY(RemoveIgTcr(TR.dds, tx2gene), tx2gene)


# count matrix for GEO

count_matrix <- counts(TR.dds.curated)

colnames(count_matrix) <- TR.dds.curated$Sample_Name

h5createFile("huHDAC7_variant/results/raw_counts.h5")

h5write(count_matrix, "huHDAC7_variant/results/raw_counts.h5", "raw_counts")

#####################################################################################################################
################################################ PCA ################################################################
#####################################################################################################################

#prepare expression values corrected for batch, to be used for PCA and heatmaps

TR.dds.curatedPCA <- rlog(TR.dds.curated, blind = F, fitType = "local")
p1 <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("Batch")) +
  ggsci::scale_color_jco()
mat <- assay(TR.dds.curatedPCA)
mat <- limma::removeBatchEffect(mat, TR.dds.curatedPCA$Batch)
assay(TR.dds.curatedPCA) <- mat
p2 <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("Treatment")) +
  ggsci::scale_color_jco()


TR.dds.curatedPCA <- vst(TR.dds.curated, blind = F, fitType = "local")
p3 <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("Batch")) +
  ggsci::scale_color_jco()
mat <- assay(TR.dds.curatedPCA)
mat <- limma::removeBatchEffect(mat, TR.dds.curatedPCA$Batch)
assay(TR.dds.curatedPCA) <- mat
p4 <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("Batch")) +
  ggsci::scale_color_jco()

pcaData <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("Batch", "Treatment"), returnData = T)

percentVar <- round(100 * attr(pcqData, "percentVar"), digits =1)

ggplot(pcaData, aes(PC1, PC2, color = Treatment, fill = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 5, shape = 22) +
  ggpubr::stat_conf_ellipse(alpha = 0.1, geom = "polygon", bary = TRUE) +
  ggpubr::stat_conf_ellipse(alpha = 1, geom = "point", bary = TRUE, level = 0.01, npoint = 2, size = 5) +
  ggsci::scale_color_jco(labels = c("control", "HDAC7 WT", "HDAC7 R166H")) +
  ggsci::scale_fill_jco()  +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_mrl() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6)) +
  guides(color = guide_legend(title= ""), fill = "none") +
  theme(legend.position = c(0.2, 0.93), aspect.ratio = 1) +
  xlim(c(-11,11)) +
  ylim(c(-11,11))
  

ggsave("huHDAC7_variant/results/figs/batchedCorrectedPCA.pdf", width = 5, height = 5)




#####################################################################################################################
#################################### Generate plot of HDAC7 expression ##############################################
#####################################################################################################################


## plot HDAC7
d <- plotCounts(TR.dds.curated, gene="ENSG00000061273", intgroup=c("Batch","Treatment"), returnData = T) %>%
  mutate(m = mean(count)) %>%
  group_by(Batch) %>%
  mutate(count = count - mean(count) + m) 

ggplot(data =d, aes(Treatment, count)) +
  stat_summary(geom = "bar", color = "black", fill = "grey") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1), width = 0.2) +
  geom_point() +
  geom_line(aes(group = Batch)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6)) +
  scale_x_discrete(labels = c("empty\nvector", "HDAC7\nWT", "HDAC7\nR166H")) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  labs(x = "")

ggsave("huHDAC7_variant/results/figs/HDAC7exp.pdf", width = 4, height = 4)


#####################################################################################################################
####################################### Run differential expression #################################################
#####################################################################################################################



degDESeqWT <- degDESeq(TR.dds.curated, Min_count = 1, alpha = 0.1, 
                          contrast = c("Treatment","hdac7WT","ctl"), fitType = "local", tx2gene = tx2gene)



degDESeqVariant <- degDESeq(TR.dds.curated, Min_count = 1, alpha = 0.1, 
                          contrast = c("Treatment","hdac7Mut","hdac7WT"), fitType = "local", tx2gene = tx2gene, lfcT = 0)



# relevel for variant vs ctl comparison
TR.dds.curated.releved <- TR.dds.curated
TR.dds.curated.releved$Treatment <- factor(TR.dds.curated.releved$Treatment, levels = c("hdac7WT", "ctl", "hdac7Mut"))


degDESeqMut2 <- degDESeq(TR.dds.curated.releved, Min_count = 1, alpha = 0.1, 
                          contrast = c("Treatment","hdac7Mut","ctl"), fitType = "local", tx2gene = tx2gene)



#####################################################################################################################
############################## Filter out high SE /low expression genes #############################################
#####################################################################################################################


# dataframe to compare lfc standard error of each gene in each comparison
lfcSEcomp <- degDESeqVariant %>% 
  mutate(Var = lfcSE) %>%
  dplyr::select(Var, gene_id, baseMean) %>%
  left_join(degDESeqMut %>% mutate(Mut = lfcSE) %>% dplyr::select(Mut, gene_id)) %>%
  left_join(degDESeqWT %>% mutate(WT = lfcSE) %>% dplyr::select(WT, gene_id))

#plot SE for each comparison
p <- ggplot(lfcSEcomp, aes(x = WT, y = Mut, color= Var)) +
  geom_point() +
  geom_abline(slope=1) +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = 2) +
  scale_colour_viridis() +
  theme(legend.position = c(0,0.8))
    
plot_grid(ggMarginal(p, fill = "grey"))

#remove genes wit SE >2 in any of the comparisons
lfcSEcomp <- lfcSEcomp %>% 
  dplyr::filter(Var < 3 & WT < 3 & Mut < 3 & baseMean >10) %>%
  dplyr::pull(gene_id)

#Keep those genes in the results tables
degDESeqVariant <- degDESeqVariant[degDESeqVariant$gene_id %in% lfcSEcomp,]
degDESeqWT <- degDESeqWT[degDESeqWT$gene_id %in% lfcSEcomp,]
degDESeqMut2 <- degDESeqMut2[degDESeqMut2$gene_id %in% lfcSEcomp,]





#####################################################################################################################
####################################### Volcano plot visualisation ##################################################
#####################################################################################################################

plot.volcano(degDESeqWT) +
  ggtitle("HDAC7 WT effect in total Tregs") +
  theme_mrl() +
  theme(legend.position = "none")

cat("\n\n")
plot.volcano(degDESeqMut2) + 
  ggtitle("HDAC7 Mut effect in total Tregs") +
  theme_mrl() +
  theme(legend.position = "none")

cat("\n\n")
plot.volcano(degDESeqVariant) + 
  ggtitle("HDAC7 variant effect in total Tregs") +
  theme_mrl() +
  theme(legend.position = "none")


#####################################################################################################################
####################################### Write DGE results to file ###################################################
#####################################################################################################################


write.table(degDESeqWT, "huHDAC7_variant/results/DGE/degDESeqWT.txt",
            quote=F,sep="\t", row.names = F, col.names = T)

write.table(degDESeqMut2, "huHDAC7_variant/results/DGE/degDESeqMut.txt",
            quote=F,sep="\t", row.names = F, col.names = T)

write.table(degDESeqVariant, "huHDAC7_variant/results/DGE/degDESeqVariant.txt",
            quote=F,sep="\t", row.names = F, col.names = T)


#####################################################################################################################
############################################# Heatmap of DEGs #######################################################
#####################################################################################################################


#Normalize counts and convert to a matrix
Norm.counts <- normTransform(TR.dds.curated)
#Norm.counts <- rlog(TR.dds.curated, blind = F, fitType = "local") 
Norm.matrix <- assay(Norm.counts)

#Select DEGs for heatmap:
degsWT <- degDESeqWT %>%
  dplyr::filter(padj < 0.1) %>%
  dplyr::select(gene_id)

degsMut <- degDESeqMut2 %>%
  dplyr::filter(padj < 0.1) %>%
  dplyr::select(gene_id)

degsVariant <- degDESeqVariant %>%
  dplyr::filter(padj < 0.1) %>%
  dplyr::select(gene_id)

#vector of gene_id for DEGs in variant comparison
degs <- degsVariant %>% 
  unique() %>%
  unname() %>%
  unlist()


#Select the DEGs in the expression matrix and replace gene_id by the gene_name
Norm.table <- as.data.frame(Norm.matrix) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(tx2gene %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%
  dplyr::filter(gene_id %in% degs) %>%
  dplyr::select(-gene_id) %>%
  tibble::column_to_rownames("gene_name")

#Convert back to matrix for heatmp input
Norm.matrix.filtered <- as.matrix(Norm.table)



#Build palette for heatmap
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(50)
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)


#replace matric column names by Treatment identity using the coldata from the DESeq object
colnames(Norm.matrix.filtered) <- as.data.frame(colData(TR.dds.curated)) %>% 
  dplyr::select(Treatment, Batch) %>% 
  unite("ID", Treatment, Batch) %>%
  unname() %>% 
  unlist()

#Split the matrix for each batch:
Norm.matrix.filtered_1 <- Norm.matrix.filtered[, 1:2]+0.1
Norm.matrix.filtered_1n <- Norm.matrix.filtered_1/rowMeans(Norm.matrix.filtered_1)

Norm.matrix.filtered_2 <- Norm.matrix.filtered[, 3:5]+0.1
Norm.matrix.filtered_2n <- Norm.matrix.filtered_2/rowMeans(Norm.matrix.filtered_2)

Norm.matrix.filtered_3 <- Norm.matrix.filtered[, 6:8]+0.1
Norm.matrix.filtered_3n <- Norm.matrix.filtered_3/rowMeans(Norm.matrix.filtered_3)

Norm.matrix.filtered_4 <- Norm.matrix.filtered[, 9:11]+0.1
Norm.matrix.filtered_4n <- Norm.matrix.filtered_4/rowMeans(Norm.matrix.filtered_4)

Norm.matrix.filtered_5 <- Norm.matrix.filtered[, 12:14]+0.1
Norm.matrix.filtered_5n <- Norm.matrix.filtered_5/rowMeans(Norm.matrix.filtered_5)

Norm.matrix.filtered_6 <- Norm.matrix.filtered[, 15:17]+0.1
Norm.matrix.filtered_6n <- Norm.matrix.filtered_6/rowMeans(Norm.matrix.filtered_6)

Norm.matrix.filtered <- cbind(Norm.matrix.filtered_1n, 
                              Norm.matrix.filtered_2n,
                              Norm.matrix.filtered_3n,
                              Norm.matrix.filtered_4n,
                              Norm.matrix.filtered_5n,
                              Norm.matrix.filtered_6n)

#Having renamed the columns by treatment we can now order it by name.
Norm.matrix.filtered <- Norm.matrix.filtered[, sort(colnames(Norm.matrix.filtered), decreasing = T)]

#Select only the WT and R166H treatment
Norm.matrix.filtered.noCTL <- Norm.matrix.filtered[,1:12]

heatmap.2(Norm.matrix.filtered, Colv=FALSE, 
          dendrogram="row", 
          labCol = colnames(Norm.matrix.filtered), 
          scale = "row", 
          col = my_palette, #cividis(50)
          symbreaks = F,
          symkey = F,
          trace = "none",
          tracecol = NA,
          lmat=rbind( c(4, 2,1), c(3,2,1), c(0,2,1)), 
          lhei=c(1, 2, 0.1), lwid = c(3,1,4),
          margins = c(11,15),
          density.info="none",
          cexRow = 0.8, cexCol = 0.8
          )


cat("\n\n")

pdf(file = "huHDAC7_variant/results/figs/heatmap.pdf", width = 7, height = 10)
heatmap.2(Norm.matrix.filtered.noCTL, Colv=FALSE, 
          dendrogram="row", 
          labCol = colnames(Norm.matrix.filtered.noCTL), 
          scale = "row", 
          col = my_palette, #cividis(50)
          symbreaks = F,
          symkey = F,
          trace = "none",
          tracecol = NA,
          lmat=rbind( c(4, 2,1), c(3,2,1), c(0,2,1)), 
          lhei=c(1, 2, 0.1), lwid = c(3,1,4),
          margins = c(11,15),
          density.info="none",
          cexRow = 0.8, cexCol = 0.8
          )
dev.off()

#####################################################################################################################
############################################# Diagonal plot #########################################################
#####################################################################################################################

# Combine data from all comparisons, renaming columns
degWT <- degDESeqWT %>%
  dplyr::select(gene_id, padj, Slfc) %>%
  plyr::rename(c("padj" = "padjWT")) %>%
  plyr::rename(c("Slfc" = "lfcWT"))

degMut <- degDESeqMut2 %>%
  dplyr::select(gene_id, padj, Slfc) %>%
  plyr::rename(c("padj" = "padjMut")) %>%
  plyr::rename(c("Slfc" = "lfcMut"))

degVariant <- degDESeqVariant %>%
  dplyr::select(gene_id, padj, Slfc, gene_name, pvalue) %>%
  plyr::rename(c("padj" = "padjVar")) %>%
  plyr::rename(c("Slfc" = "lfcVar")) %>%
  plyr::rename(c("pvalue" = "pvalueVar"))

CombinedDEGs <- full_join(degWT, degMut, by = "gene_id")
CombinedDEGs <- full_join(CombinedDEGs, degVariant, by = "gene_id")

CombinedDEGsplots <- CombinedDEGs %>% 
  mutate(FDR_cat = ifelse(padjVar < 0.1, "FDR < 0.1", "FDR > 0.1")) %>%
  mutate(lfcWT = ifelse(lfcWT > 4, 4, ifelse(lfcWT < -4, -4, lfcWT))) %>%
  mutate(lfcMut = ifelse(lfcMut > 4, 4, ifelse(lfcMut < -4, -4, lfcMut))) %>%
  mutate(gendegs = ifelse(padjWT<0.1, "Differentially expressed vs control", ifelse(padjMut<0.1, "Differentially expressed vs control", "not differentially expressed"))) %>%
  unite("combined",gendegs, FDR_cat, sep = ".", remove = FALSE) %>%
  mutate(FDR_cat = ifelse(combined == "Differentially expressed vs control.FDR > 0.1", "control FDR < 0.1", FDR_cat))


# build plot
diffplot <- ggplot(CombinedDEGsplots, aes(x = lfcWT, y = lfcMut)) + # sets general plot structure
  geom_point(fill = "dark grey", alpha = 1, na.rm = T, size = 2, shape = 21, color = "dark grey") +
  stat_density2d(aes(fill = stat(level)), bins = 100, geom = "polygon") +
  geom_point(data = subset(CombinedDEGsplots, FDR_cat %in% c("FDR < 0.1")), size = 3, fill = "white", alpha = 1, na.rm = T, shape = 21, stroke = 0) +
  geom_point(data = subset(CombinedDEGsplots, FDR_cat %in% c("FDR < 0.1") & lfcVar > 0) , size = 3, alpha = 1, na.rm = T, fill = "red4", shape = 21, stroke = 0.2) +
  geom_point(data = subset(CombinedDEGsplots, FDR_cat %in% c("FDR < 0.1") & lfcVar < 0) , size = 3, alpha = 1, na.rm = T, fill = "navy", shape = 21, stroke = 0.2) +
  #geom_point(data = subset(CombinedDEGsplots, FDR_cat %in% c("FDR < 0.1")) , size = 3, alpha = 1, na.rm = T, fill = NA, shape = 21, stroke = 0.2) +
  theme_minimal(base_size = 16) + # sets visual theme of the plot (background, fonts etc.)
  theme(legend.position = "none") + # legend formatting
  #geom_smooth( method = "loess", colour = "#90d092", alpha = 0.5, fill = "#90d092") +
  theme(aspect.ratio = 1) +
  labs(x = "LFC WT vs control", y = "LFC R166H vs control") +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = 0.5, linetype = "dotted") +
  geom_abline(slope = 1, intercept = -0.5, linetype = "dotted") +
  guides( color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  scale_fill_gradient(low = "dark grey", high = "black")

### Select 100 genes with lowest FDR:
top_labelled <- top_n(CombinedDEGsplots, n = 100, wt = (padjVar*-1)) %>%
  dplyr::filter(grepl("RSAD2", gene_name) |
                  grepl("HIF1A", gene_name) |
                  grepl("TCF4", gene_name) |
                  grepl("ATRIP", gene_name) |
                  grepl("LRRC29", gene_name) |
                  grepl("PIP", gene_name) |
                  grepl("SLAMF7", gene_name) |
                  grepl("CD226", gene_name))


## Labels dots corresponiding to 10 lowest FDR genes with gene names:
diffplot<-diffplot + geom_label_repel(data = top_labelled,
                                      mapping = aes(label = gene_name, color = lfcVar < 0),
                                      size = 3.5,
                                      label.size = 0.2,
                                      segment.size = 0.2,
                                      #force = 300,
                                      #ylim = c(-1,1.5,1.5),
                                      #xlim = c(-1.5,1.5,1),
                                      nudge_x= 1.3*top_labelled$lfcWT,
                                      nudge_y= 1.3*top_labelled$lfcMut,
                                      max.iter = 5000,
                                      fontface = 'bold',
                                      #color = 'black',
                                      box.padding = unit(0.5, "lines"),
                                      point.padding = unit(0.5, "lines")) +
  scale_color_manual(values = c("red4", "navy"))

# fine tune appearance
diffplot + 
  xlim(-1.5,1.5) + 
  ylim(-1.5,1.5) + 
  theme_mrl(1) +
  theme(axis.ticks = element_line(),
        axis.line = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        aspect.ratio = 1,
        legend.key.size = unit(0.3, "cm"))




ggsave("huHDAC7_variant/results/figs/diagonal_plot_2.pdf",
       width = 5,
       height = 5)


