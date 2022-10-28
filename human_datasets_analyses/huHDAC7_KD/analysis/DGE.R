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
all.samples <- read_csv("huHDAC7_KD/data/GEO_table_STM_huKD.csv") %>% spread(fastq_read, file)

all.samples$Batch <- as.factor(all.samples$Batch)
all.samples$Treatment <- factor(all.samples$Treatment, levels = c("shHDAC7", "scrM"))
all.samples$sampleID <- as.factor(all.samples$sampleID)

#modify sample table to include quantitations files paths for each samples (gtf and results files)
sample.table.RNA <- all.samples %>%
  dplyr::filter(!is.na(RNA_FastQ1) & !is.na(RNA_FastQ2)) %>%
  mutate(RNA_Htseq_File = paste0("huHDAC7_KD/results/rna/quantitation/htseq-count/", Sample_Name, ".gtf"),
         RNA_RSEM_File = paste0("huHDAC7_KD/results/rna/quantitation/rsem/", Sample_Name, ".genes.results"))

# Import RSEM-quantitated reads with tximport:
txi.rsem <- tximport(sample.table.RNA$RNA_RSEM_File,
                     type = "rsem", countsFromAbundance = "lengthScaledTPM")
#prevents error when importing in DESeq
txi.rsem$length[txi.rsem$length == 0] <- 1

#Import data from rsem quantitation into a dds object
main.dds <- DESeqDataSetFromTximport(txi.rsem, all.samples, design = ~1)


#Change the design after removing the deprecated levels of cell type and treatment factors.
#Mandatory to input the design only at the end otherwise the design matrix will be of the wrong size.
main.dds@design <- ~sampleID+Treatment


#Remove irrelevant gene entries: genes on the Y chromosome and Ig/TCR segments
main.dds.curated.M <- RemoveY(RemoveIgTcr(main.dds, tx2gene), tx2gene)



# count matrix for GEO

library(rhdf5)
count_matrix <- counts(main.dds.curated.M)

colnames(count_matrix) <- main.dds.curated.M$Sample_Name

h5createFile("huHDAC7_KD/results/raw_counts.h5")

h5write(count_matrix, "huHDAC7_KD/results/raw_counts.h5", "raw_counts")

#####################################################################################################################
################################################ PCA ################################################################
#####################################################################################################################

#plot PCA corrected per donor

TR.dds.curatedPCA <- vst(main.dds.curated.M, blind = F, fitType = "local")
p3 <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("sampleID")) +
  ggsci::scale_color_jco() +
  theme_mrl(0.5) +
  theme(legend.position = "bottom", aspect.ratio = 1)
mat <- assay(TR.dds.curatedPCA)
mat <- limma::removeBatchEffect(mat, TR.dds.curatedPCA$sampleID)
assay(TR.dds.curatedPCA) <- mat


## Plot after correction
p4 <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("Treatment")) +
  ggsci::scale_color_jco() +
  theme_mrl(0.5) +
  theme(legend.position = "bottom", aspect.ratio = 1)
plot_grid(p3,p4)
pcqData <- DESeq2::plotPCA(TR.dds.curatedPCA, intgroup = c("sampleID", "Treatment"), returnData = T)

percentVar <- round(100 * attr(pcqData, "percentVar"), digits =1)

ggplot(pcqData, aes(PC1, PC2, color = Treatment, fill = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 5, shape = 22) +
  ggpubr::stat_conf_ellipse(alpha = 0.1, geom = "polygon", bary = TRUE) +
  ggpubr::stat_conf_ellipse(alpha = 1, geom = "point", bary = TRUE, level = 0.01, npoint = 2, size = 5) +
  ggsci::scale_color_jco() +
  ggsci::scale_fill_jco()  +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_mrl() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6)) +
  guides(color = guide_legend(title= ""), fill = "none") +
  theme(legend.position = c(0.2, 0.93), aspect.ratio = 1) +
  xlim(c(-11,11)) +
  ylim(c(-11,11))


ggsave("huHDAC7_KD/results/figs/batchedCorrectedPCA.pdf", width = 5, height = 5)




#####################################################################################################################
#################################### Generate plot of HDAC7 expression ##############################################
#####################################################################################################################


## plot HDAC7
d <- plotCounts(main.dds.curated.M, gene="ENSG00000061273", intgroup=c("sampleID","Treatment"), returnData = T) %>%
  mutate(Treatment = factor(Treatment, levels = c("scrM", "shHDAC7")))

ggplot() +
  stat_summary(data =d, aes(Treatment, count), geom = "bar") +
  geom_point(data =d, aes(Treatment, count)) +
  geom_line(data =d, aes(Treatment, count, group = sampleID)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6)) +
  scale_x_discrete(labels = c("Scrambled\nshRNA", "HDAC7\nshRNA")) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
  labs(x = "")

ggsave("huHDAC7_KD/results/figs/HDAC7exp.pdf", width = 3, height = 4)



#####################################################################################################################
####################################### Run differential expression #################################################
#####################################################################################################################



# DGE analysis
main.dds.curated.M@design <- ~sampleID+Treatment

degDESeqH7m <- degDESeq(main.dds.curated.M, Min_count = 1, alpha = 0.1, 
                        contrast = c("Treatment","shHDAC7", "scrM"), fitType = "local", tx2gene = tx2gene)

degDESeqH7m <- degDESeqH7m %>% dplyr::filter(lfcSE < 2)


p1 <- plot.volcano(degDESeqH7m, top = 70) + 
  scale_y_continuous(limits = c(15,61)) +
  labs(title = "Effect of HDAC7 knock-down on gene expression in Tregs") +
  theme_mrl(0.5) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none", axis.line.x = element_blank(), axis.title.y = element_blank())

p2 <- plot.volcano(degDESeqH7m, top = 80) + 
  scale_y_continuous(limits = c(0,15)) +
  labs(caption = paste0(nrow(degDESeqH7m %>% dplyr::filter(padj < 0.1 & Slfc < 0)), 
                        " genes downregulated and ",
                        nrow(degDESeqH7m %>% dplyr::filter(padj < 0.1 & Slfc > 0)),
                        " genes upregulated")) +
  theme_mrl(0.5) +
  theme(plot.margin = margin(t=0, unit="cm"),
        legend.position = c(0.9,0.3)) +
  guides(color=guide_legend(title="Significance"))

cat("\n\n")

plot_grid(p1, p2, ncol = 1, axis = "lr", align = "v", rel_heights = c(2,3))





colnames(degDESeqH7m) <- paste(colnames(degDESeqH7m), "H7m", sep = "_")


### plot MA plot


genes <- degDESeqH7m %>% top_n(n=30, wt = abs(Slfc_H7m)) %>% pull(gene_name_H7m)
plotH7m <- annotated.plot(degDESeqH7m %>% 
                            mutate(X = log2(baseMean_H7m)) %>%
                            mutate(Y = Slfc_H7m) %>%
                            mutate(FDR_cat = ifelse(padj_H7m >= 0.1 | is.na(padj_H7m), 
                                                    "FDR > 0.1", "FDR < 0.1")) %>%
                            mutate(gene_name = gene_name_H7m),
                          genes,
                          bins = 8, nudgeX = 0, nudgeY = 0.2) 



plotH7m[[2]] +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 0),
        axis.ticks.x = element_line(color = "black", size = 1/2, linetype = "solid"),
        axis.ticks.y = element_line(color = "black", size = 1/2, linetype = "solid")) +
  labs(x = "log10 mean of normalized counts", y = "log2 fold change") +
  lims(y = c(-3.5,3.5))

ggsave("huHDAC7_KD/results/figs/MAplot.pdf", width = 7, height = 4)


#####################################################################################################################
############################################# Heatmap of DEGs #######################################################
#####################################################################################################################


#Normalize counts and convert to a matrix
Norm.counts <- normTransform(main.dds.curated.M)
Norm.matrix <- assay(Norm.counts)

#Select DEGs for heatmap:
degs <- degDESeqH7m %>%
  dplyr::filter(padj_H7m < 0.1) %>%
  dplyr::select(gene_id_H7m) %>%
  unique() %>%
  unname() %>%
  unlist()


#Select the DEGs in the expression matrix and replace gene_id by the gene_name
Norm.table <- as.data.frame(Norm.matrix) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(tx2gene %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%
  dplyr::filter(gene_id %in% degs) %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
  dplyr::select(-gene_id) %>%
  tibble::column_to_rownames("gene_name")


rowLabels <- data.frame(names = as.character(rownames(Norm.table))) %>%
  dplyr::mutate(names = ifelse(grepl("IFN|HIF|CXCR|CCR|CCL|CXCL|^IL|TIGIT|CD226|FOXO|TBX21|FOXP3|IL10|IFNG|LAG3|HAVCR2|PDCD1|IFNGR|SP1R|PIP4K|TCF4|CTNNB1|CTLA4|IL2R|S1PR1|IL21R|CD74|CD47|CD5|ITGB7|ITGA4|IFITM|IRF|BTLA|IL2RG|CD24|CXCL8|IL26|CEBPA|IGF1|IL1B|CCR2|ITGA1", names, ignore.case = T), as.character(names), ""))
#Convert back to matrix for heatmp input
Norm.matrix.filtered <- as.matrix(Norm.table)

# Venn Diagram of 3 comparisons overlap
#Venn <- list(CTLvsWT = degsWT %>% unname() %>% unlist(),
#            CTLvsR166H = degsMut %>% unname() %>% unlist(),
#           WTvsR166H = degsVariant %>% unname() %>% unlist())

#Venn.plot <- venn.diagram(Venn, "Venn.tiff")

#Build palette for heatmap
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(50)
my_palette <- colorRampPalette(c("royalblue3", "black", "gold"))(n = 299)


#replace matric column names by Treatment identity using the coldata from the DESeq object
colnames(Norm.matrix.filtered) <- as.data.frame(colData(main.dds.curated.M)) %>% 
  dplyr::select(Treatment, sampleID) %>% 
  unite("ID", Treatment, sampleID) %>%
  unname() %>% 
  unlist()

#Having renamed the columns by treatment we can now order it by name.
Norm.matrix.filtered <- Norm.matrix.filtered[, sort(colnames(Norm.matrix.filtered), decreasing = F)]

# split the matrix to normalize per donor

Norm.matrix.filtered.41 <- Norm.matrix.filtered[,grepl("1", colnames(Norm.matrix.filtered))]+0.1
Norm.matrix.filtered.42 <- Norm.matrix.filtered[,grepl("2", colnames(Norm.matrix.filtered))]+0.1
Norm.matrix.filtered.43 <- Norm.matrix.filtered[,grepl("3", colnames(Norm.matrix.filtered))]+0.1

Norm.matrix.filtered.41 <- Norm.matrix.filtered.41/rowMeans(Norm.matrix.filtered.41)
Norm.matrix.filtered.42 <- Norm.matrix.filtered.42/rowMeans(Norm.matrix.filtered.42)
Norm.matrix.filtered.43 <- Norm.matrix.filtered.43/rowMeans(Norm.matrix.filtered.43)

# merge normalized matrix back together
Norm.matrix.filtered <- cbind(Norm.matrix.filtered.41,
                              Norm.matrix.filtered.42,
                              Norm.matrix.filtered.43)



# plot heatmap

# Uncomment to save as pdf
pdf(file="huHDAC7_KD/results/figs/Heatmap_STM.pdf", width = 4, height = 10)

#plotting functon
hm <- heatmap.2(Norm.matrix.filtered, Colv=T,
                labRow = rowLabels$names,
                dendrogram="row", 
                labCol = colnames(Norm.matrix.filtered), 
                hclustfun = function(x) hclust(x,method = "average"),
                scale = "row", 
                col = my_palette, #cividis(50)
                #ColSideColors = c("grey","grey","grey","black","black","black"),
                #RowSideColors = Order,
                symbreaks = F,
                symkey = F,
                trace = "none",
                tracecol = NA,
                # lmat=rbind( c(0,0,0,0,0),
                #             c(0,0,5,5,0),
                #             c(0,6,6,2,0),
                #             c(0,4,1,3,0),
                #             c(0,0,0,0,0)),
                # lhei=c(0.1,0.5,1,4,0.1), lwid = c(0.1,1,0.2,4,0.1),
                # lmat=rbind( c(0,0,0,0,0),
                #             c(4,4,2,1,0),
                #             c(0,3,2,1,0),
                #             c(0,3,2,1,0),
                #             c(0,3,2,1,0),
                #             c(0,0,2,1,0),
                #             c(0,0,0,0,0)),
                # lhei=c(0.1,2, 2, 2, 2, 0.1,1), lwid = c(1,2,2,4,1),
                # margins = c(1,1),
                density.info="none",
                cexRow = 0.5, cexCol = 0.8)
# Uncomment to save as pdf
dev.off()



write_csv(degDESeqH7m, "huHDAC7_KD/results/DGE/huHDAC7_KD_STM.csv")
