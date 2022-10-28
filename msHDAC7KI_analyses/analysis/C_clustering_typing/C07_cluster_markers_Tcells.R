## For R3.5!


#################################################
################### Set up ######################
#################################################

## Take the input argument
# args <- commandArgs(TRUE)
# sample_name <- args[1]
# exp_matrix <- args[2]
# demux_best <- args[3]
# out_path <- args[4]

## load libraries

library(gridExtra)
library(knitr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggExtra)
library(stats)
# library(RColorBrewer)
library(viridis)
library(scales)
library(ggrepel)
library(circlize)
library(RColorBrewer)
library(gplots)
library(cowplot)
library(gtable)
library(data.table)
library(Seurat)
library(hdf5r)
library(ggsci)
library(harmony)
library(presto)



## ggplot theme 

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






#########################################################
################### paths and stems #####################
#########################################################



# Set wd on project folder
setwd("/home/pa326/project/HDAC7_KI_EAE_10x/")


results_dir <- "results/master"
dir.create(results_dir, showWarnings = FALSE)


objects_dir <- file.path(results_dir,"seurat_objects")
dir.create(objects_dir, showWarnings = FALSE)



results_dir <- file.path(results_dir,"clustering_typing")
dir.create(results_dir, showWarnings = FALSE)

#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- readRDS(file.path(objects_dir,"master_seurat_harmony_Tcells_low_res_clustered.rds"))

# manually select appropriate sublclustering resolutions
meta_clustering <- merged.seurat@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  mutate(SNNk20_algo1_res0.3 = as.character(SNNk20_algo1_res0.3)) %>%
  mutate(Sub_res0.15 = as.character(Sub_res0.15)) %>%
  mutate(Sub_res0.1 = as.character(Sub_res0.1)) %>%
  mutate(curated_cluster = ifelse(SNNk20_algo1_res0.3 %in% c("0","1"), Sub_res0.15, Sub_res0.1)) %>%
  tibble::column_to_rownames("cell.id")

merged.seurat@meta.data <- meta_clustering

cat("breakdown of curated clusters")
table(meta_clustering$curated_cluster)

############################
####### find markers ####### 
############################


#----------------------------
# general DGE for each cluster


presto_results <- wilcoxauc(merged.seurat, "curated_cluster", assay = 'data')

saveRDS(presto_results, file.path(results_dir,"presto_markers_Tcells_curated_cluster.rds"))


#----------------------------
# DGE paired for each cluster


cluster_assignments <- merged.seurat@meta.data$curated_cluster

combinations <- combn(unique(cluster_assignments),
                        2,
                        simplify = F)
                        
paired_results <- lapply(combinations, function(comb) {
    wilcoxauc(merged.seurat, "curated_cluster", groups_use = comb, assay = 'data') %>%
        mutate(groupA = comb[1]) %>%
        mutate(groupB = comb[2])
})

saveRDS(paired_results, file.path(results_dir,"presto_paired_markers_Tcells_curated_cluster.rds"))



paired_results <- do.call(rbind, paired_results) %>%
  dplyr::filter(padj < 0.01 & logFC > 0 & auc > 0.6)
  
saveRDS(paired_results, file.path(results_dir,"presto_paired_markers_Tcells_curated_cluster_filtered.rds"))


