#####################################################################################################################
############################################## Set up ###############################################################
#####################################################################################################################


##############
## packages ##
##############
library(gridExtra)
library(knitr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(stats)
# library(RColorBrewer)
library(viridis)
library(scales)
library(ggrepel)
library(EnsDb.Hsapiens.v86)
#library(DEFormats)
#library(statmod)
#library(Cairo)
#library(vsn)
#library(kableExtra)
library(gplots)
library(cowplot)
library(ggdendro)
library(gtable)
library(Seurat)
library(ggsignif)
library(forcats)
library(lmerTest)


knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
#knitr::opts_chunk$set(dev = 'png', dpi = 600)

#set root directory to be one level up of the current directory
#knitr::opts_knit$set(root.dir = '..')


WRITE_FILES <- FALSE


##################
## ggplot theme ############################################################################
##################

theme_mrl <- function(x = 1) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(size = 12*x,face = "bold", angle = 0, vjust = 0.6),
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




#####################################################################################################################
############################################## Seurat pipeline ######################################################
#####################################################################################################################

#################################
## data retrieval and clean up ##
#################################

#Datasets from Miragaia et al. Immunity 2019
load("sr_list_SS2.RData")
human_seurat <- sr_list$human

#Set of co-expressed gene from PPA
HDAC7_module <- read.table("HDAC7moduleTreg_wID.txt", header = TRUE, sep = " ", dec = ".", fill = TRUE)

#Gene-Ensembl IDs annotations file
edb <- EnsDb.Hsapiens.v86
tx2gene <- genes(edb, return.type = "data.frame")

#We export the slots with raw counts, normalized data and scaled data from the Seurat object and attempt to create respective assays slots.
raw_counts <- as.matrix(human_seurat@raw.data)
dim(raw_counts)
head(rownames(raw_counts))
head(colnames(raw_counts))

norm_data <- as.matrix(human_seurat@data)
dim(norm_data)
head(rownames(norm_data))
head(colnames(norm_data))

scaled_data <- as.matrix(human_seurat@scale.data)
dim(scaled_data)
head(rownames(scaled_data))
head(colnames(scaled_data))

#We also export the metadata table
metadata <- as.data.frame(human_seurat@meta.data)

#We create a new seurat object initialized with the raw counts
new_seurat <- CreateSeuratObject(counts = raw_counts, meta.data = metadata)

#We add in the normalized and scaled data
new_seurat@assays$RNA@data <- norm_data
new_seurat@assays$RNA@scale.data <- scaled_data



###################################
## dim reduction and integration ##
###################################

Idents(new_seurat) <- "tissue_cell"


#We identify variable features, meaning features that are outliers on a 'mean variability plot'.
new_seurat <- FindVariableFeatures(object = new_seurat, selection.method = "mean.var.plot", nfeatures = 2000)
new_seurat <- RunPCA(object = new_seurat, npcs = 50, verbose = TRUE)

#We integrate the datasets by tissue, as there seems to be a prominent tissue effect on a basic parameter such as the number of unique genes.
blood <- subset(x = new_seurat, subset = tissue == "blood")
skin <- subset(x = new_seurat, subset = tissue == "skin")


tissue.anchors <- FindIntegrationAnchors(object.list = list(blood, skin), dims = 1:20)
seurat.combined <- IntegrateData(anchorset = tissue.anchors, dims = 1:20)

#scaling and PCA
DefaultAssay(object = seurat.combined) <- "integrated"

seurat.combined <- ScaleData(object = seurat.combined, verbose = TRUE)

seurat.combined <- RunPCA(object = seurat.combined, npcs = 50, verbose = TRUE)

#We select the number of significant PCs to be used in the downstream clustering based on an elbow plot
ElbowPlot(seurat.combined)



#############################################################
## Clustering steps for tSNE, UMAP and neighrest neighbour ##
#############################################################

seurat.combined <- subset(seurat.combined, subset = cell.type == "Treg")
seurat.combined <- RunTSNE(object = seurat.combined, reduction = "pca", dims = 1:10)
seurat.combined <- RunUMAP(object = seurat.combined, reduction = "pca", dims = 1:10)
seurat.combined <- FindNeighbors(object = seurat.combined, reduction = "pca", k.param = 20, dims = 1:10, verbose = TRUE)
seurat.combined <- FindClusters(object = seurat.combined, resolution = 0.7)


##########################
## HDAC7 module scoring ##
##########################

#Set of co-expressed genes
HDAC7_module <- read.table("huHDAC7_coexpression/results/HDAC7moduleTreg_wID.txt", header = TRUE, sep = " ", dec = ".", fill = TRUE)


HDAC7 <- list(c(as.character(HDAC7_module$gene_id)))
seurat.combined <-AddModuleScore(object = seurat.combined, features = HDAC7, assay = "RNA", name = "HDAC7_coexpression_module")

saveRDS(seurat.combined, "huHDAC7_coexpression/results/TeichTreg/Teichman.Treg.skinblood.seurat.rds")

# tissue effect test with mixed model to account for cells nested in donor
lmerTest::lmer(HDAC7_coexpression_module1 ~ tissue + (1 | individual), seurat.combined.reclust@meta.data) %>% summary()


#####################################################################################################################
######################################## Plot results and save ######################################################
#####################################################################################################################

########################
## add cluster labels ##
########################


# add curated labels to clusters
new.meta <- seurat.combined.reclust@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  mutate(cluster_labels = ifelse(integrated_snn_res.0.7 == "0", "2_tissue-resident",
                                 ifelse(integrated_snn_res.0.7 == "1", "0_circulating",
                                        ifelse(integrated_snn_res.0.7 == "2", "1_intermediate", "3_stressed"))))


# add back to seurat object
seurat.combined.reclust@meta.data <- new.meta %>%
  tibble::column_to_rownames("cell.id")

################
## UMAP plots ##
################

# Cluster by color
DimPlot(object = seurat.combined.reclust, reduction = "umap", pt.size = 5, group.by = "cluster_labels") + theme(aspect.ratio = 1.2) +
  ggsci::scale_color_futurama() +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", fill=NA, size= 0.5), 
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #legend.text = element_text(size = 5),
        legend.position = "right",
        legend.box.margin =margin(),
        #legend.justification = c(0.3,0.8),
        aspect.ratio = 1)

ggsave("huHDAC7_coexpression/results/figs/clusterUMAP_STM.pdf", width = 5, height = 5)

DimPlot(object = seurat.combined.reclust, reduction = "umap", pt.size = 2, group.by = "tissue") + theme(aspect.ratio = 1.2) +
  scale_color_manual(values = ggsci::pal_jco()(10)[c(4,6)]) +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", fill=NA, size= 0.5), 
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 5),
        legend.position = c(0,1),
        legend.box.margin =margin(),
        legend.justification = c(0.25,0.8),
        aspect.ratio = 1)

ggsave("huHDAC7_coexpression/results/figs/tissueUMAP_STM.pdf", width = 2, height = 2)


#################
## module plot ##
#################



p2 <-VlnPlot(seurat.combined, reduction = "tsne",  features = "HDAC7_coexpression_module1", coord.fixed = TRUE, group.by = "tissue", pt.size = 0 ) +
  geom_signif(comparisons = list(c("blood", "skin")),
              map_signif_level = TRUE, textsize=6) +
  ylim(0, 0.2) +
  scale_fill_manual(values = ggsci::pal_jco()(10)[c(4,6)]) +
  #geom_boxplot() +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", fill=NA, size= 0.5), 
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none",
        aspect.ratio = 2) +
  labs(x = "", y = "Module score", title = "Module n 1\nsignature") +
  geom_boxplot()


p2$layers[[1]] <- NULL
p2


ggsave("huHDAC7_coexpression/results/figs/HDAC7moduleScore_STM.pdf", width = 2, height = 5)



################################
#### split between tissues #####
################################

seurat.combined.reclust@meta.data %>%
  ggplot(aes(fct_rev(cluster_labels), fill = tissue)) +
  geom_bar() +
  scale_y_continuous(expand = expand_scale(add = c(0,2))) +
  scale_fill_manual(values = ggsci::pal_jco()(10)[c(4,6)]) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        aspect.ratio = 1,
        axis.title = element_blank()) +
  coord_flip()

ggsave("huHDAC7_coexpression/results/figs/cluster_tissue_plot.pdf", width = 5, height = 5)


seurat.combined.reclust@meta.data %>%
  ggplot(aes(fct_rev(cluster_labels), fill = tissue)) +
  geom_bar(position = "fill") +
  scale_y_continuous(expand = expand_scale(add = c(0,2))) +
  scale_fill_manual(values = ggsci::pal_jco()(10)[c(4,6)]) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        aspect.ratio = 1,
        axis.title = element_blank()) +
  coord_flip()


p2 <-VlnPlot(seurat.combined.reclust, reduction = "tsne",  features = "HDAC7_coexpression_module1", coord.fixed = TRUE, group.by = "cluster_labels", pt.size = 0) +
  geom_signif(comparisons = list(c("2_tissue-resident", "1_intermediate"), c("0_circulating","1_intermediate"),
                                 c("0_circulating","2_tissue-resident")),
              step_increase = 0.03,
              map_signif_level = F, textsize=2) +
  ylim(0, 0.2) +
  #scale_fill_manual(values = ggsci::pal_jco()(10)[c(4,6)]) +
  #geom_boxplot() +
  theme(panel.border = element_rect(colour = "black", linetype = "solid", fill=NA, size= 0.5), 
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        aspect.ratio = 2) +
  labs(x = "", y = "Module score", title = "Module n 1\nsignature") +
  geom_boxplot()  +
  ggsci::scale_fill_futurama()


p2$layers[[1]] <- NULL
p2


ggsave("huHDAC7_coexpression/results/figs/HDAC7moduleScore_cluster_STM.pdf", width = 5, height = 5)

#####################
## dot plot seurat ##
#####################

gene_ids <- tx2gene %>%
  dplyr::filter(seq_name %in% c(1:22, "X", "Y")) %>%
  dplyr::filter(gene_name %in% c("ARF6", "IQGAP2", "TAOK3", "RORA", "CTLA4", "CREM", "LTB"))


dotplot <-
  ggplot_build(DotPlot(seurat.combined.reclust, gene_ids$gene_id, assay = "RNA", group.by = "cluster_labels"))$plot$data %>%
  dplyr::rename("gene_id" = "features.plot") %>%
  dplyr::filter(id != "3") %>%
  left_join(gene_ids %>% dplyr::select(gene_id, gene_name)) %>%
  mutate(gene_name = factor(gene_name, levels = c("CREM", "CTLA4", "RORA", "ARF6", "IQGAP2", "LTB", "TAOK3")))

ggplot(dotplot, aes(gene_name,fct_rev(id), size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21) +
  scale_fill_viridis() +
  #scale_fill_viridis(limits = quantile(dotplot$avg.exp.scaled, c(0.025,0.975)), oob = squish) +
  scale_size(range = c(1,7)) +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", linetype = "solid", fill=NA, size= 0.5),
        legend.text =  element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.box ="horizontal",
        aspect.ratio = 1)


ggsave("huHDAC7_coexpression/results/figs/dotplot.pdf", width = 5, height = 5)
