#################################################
################### Set up ######################
#################################################









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
library(ggridges)
library(sp)
library(harmony)


theme_mrl <- function(x = 1) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(size = 12*x,face = "bold", angle = 45, vjust = 0.9, hjust=0.9),
      axis.text.y = element_text(size = 12*x,face = "bold"),
      axis.title.x = element_text(size = 12*x,face = "bold"),
      axis.title.y = element_text(size = 12*x,face = "bold"),
      strip.background = element_rect(fill="gray20", colour="gray20", linetype="solid"),
                                      strip.text = element_text(size=12*x, colour="white", face="bold"),
                                      legend.title = element_text(size=14*x, face = "bold"),
                                      legend.text = element_text(size=12*x, color="gray20", face="bold"),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
      plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
      plot.caption = element_text(hjust = 0, face = "italic")
                                      )
}





#################################################
################ Import data ####################
#################################################


# Set wd on project folder
setwd("/home/pa326/project/HDAC7_KI_EAE_10x/")


results_dir <- "results/master"
dir.create(results_dir, showWarnings = FALSE)


objects_dir <- file.path(results_dir,"seurat_objects")
dir.create(objects_dir, showWarnings = FALSE)


# import seurat object

merged.seurat <- readRDS(file.path(objects_dir,"master_seurat_QCed.rds"))


# filter out unwanted cells from first pass analysis: removing myeloid, B cells, damaged/stress T cells


keep.cells <- readRDS(file.path(results_dir, "clustering_typing", "master_seurat_low_res_clustering_table_harmony.rds")) %>%
    dplyr::filter(!grepl("2_1|^[4-6]",Sub_res0.05)) %>%
    pull(cell.id)
    
    
keep.cells <- rownames(merged.seurat@meta.data) %in% keep.cells


cat("\n breakdown of cells to keep/remove \n")
summary(keep.cells)


merged.seurat <- merged.seurat[,keep.cells]


#########################################################
################# Input/output info #####################
#########################################################


results_dir <- file.path(results_dir,"reductions")
dir.create(results_dir, showWarnings = FALSE)


#output name stem
out_stem <- "harmony_Tcells"



#########################################################
############ Count normalization & PCA ##################
#########################################################



#normalize for count/cell and log transform
merged.seurat <- NormalizeData(merged.seurat)

#Variable feature selection for PCA
merged.seurat <- FindVariableFeatures(merged.seurat, selection.method = "vst", nfeatures = 2000)
p1 <- VariableFeaturePlot(merged.seurat) + 
    labs(caption = "hypervariable genes selection for PCA (defaults to 2000)") + 
    theme(legend.position = "none") +
    geom_density2d() +
    scale_y_log10()

# for debugging, sample a few cells from the main object
#merged.seurat <- merged.seurat[, sample(1:nrow(init.merged.seurat@meta.data),1000)]



# Filter out unwanted genes.


HVG <- merged.seurat@assays[["RNA"]]@meta.features %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("gene_name") %>%
                    dplyr::filter(vst.variance.standardized >1) %>% # removed  & vst.mean > 10^-3
                    dplyr::filter(gene_category == "other" | grepl("C_gene$", gene_biotype)) %>% # keep constant chain from natigen receptor gene, might help for distinguishing gd and ab T cells
                    pull(gene_name)

cat("\nnumber of HVG", length(HVG))

# Run and plot PCA
merged.seurat <- ScaleData(merged.seurat, features = HVG)
merged.seurat <- RunPCA(merged.seurat, features = HVG)



cat("\nmerging meta and PCA embeddings\n")
reductions.meta <- merged.seurat@reductions$pca@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
		left_join(merged.seurat@meta.data %>% tibble::rownames_to_column("cell.id"), by = "cell.id") 

cat("\ndim(pca.meta)\n")
dim(reductions.meta)


p2 <- ElbowPlot(merged.seurat, ndims =50) + labs(title = paste("PCA on", length(HVG), "HVGs"),
                                                caption = "Variance explained by each PC")

#Tracy Widom test
tw.stat <- merged.seurat@reductions$pca@stdev
tw.stat <- AssocTests::tw(tw.stat^2, length(tw.stat))$statistic

p3 <- as.data.frame(tw.stat) %>%
  tibble::rownames_to_column("TW") %>%
  mutate(PC = as.numeric(sub("TW","",TW))) %>%
  ggplot(aes(PC, tw.stat)) +
  geom_point() +
  geom_hline(yintercept = 3.2724, linetype = "dotted") +
  geom_hline(yintercept = 0.9793, linetype = "dotted") +
  geom_text(aes(50,3.2724,label = "p=0.001", vjust = -0.5, hjust = 1)) +
  geom_text(aes(50,0.9793,label = "p=0.05", vjust = -0.5, hjust = 1)) +
  geom_text(aes(label = PC, vjust = -1), check_overlap = T) +
  labs(y = "Tracy Widom statistic", caption = "Tracy Widom test")


# PC loadings
p5 <- VizDimLoadings(merged.seurat, dims = c(1), reduction = "pca") + labs(caption = "loadings for first PC")
p6 <- VizDimLoadings(merged.seurat, dims = c(2), reduction = "pca") + labs(caption = "loadings for 2nd PC")
p7 <- VizDimLoadings(merged.seurat, dims = c(3), reduction = "pca") + labs(caption = "loadings for 3rd PC")
p8 <- VizDimLoadings(merged.seurat, dims = c(4), reduction = "pca") + labs(caption = "loadings for 4th PC")



ggsave(file.path(results_dir,"master_seurat_PCA_Tcells.pdf"), plot_grid(p1,p2,p3,p5,p6,p7,p8,ncol = 3), width = 16, height = 16)

saveRDS(reductions.meta, file.path(results_dir,"master_seurat_meta_table_Tcells.rds"))


#########################################################
################### harmony 1st pass ####################
#########################################################

merged.seurat <- merged.seurat %>% RunHarmony("sample", plot_convergence = TRUE) 


print("merging meta and PCA embeddings")
pca.meta <- merged.seurat@reductions$harmony@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
		left_join(merged.seurat@meta.data %>% tibble::rownames_to_column("cell.id"), by = "cell.id") 




p2 <- ElbowPlot(merged.seurat, ndims =50, reduction = "harmony") + labs(caption = "Variance explained by each PC")

#Tracy Widom test
tw.stat <- merged.seurat@reductions$harmony@stdev
tw.stat <- AssocTests::tw(tw.stat^2, length(tw.stat))$statistic

p3 <- as.data.frame(tw.stat) %>%
  tibble::rownames_to_column("TW") %>%
  mutate(PC = as.numeric(sub("TW","",TW))) %>%
  ggplot(aes(PC, tw.stat)) +
  geom_point() +
  geom_hline(yintercept = 3.2724, linetype = "dotted") +
  geom_hline(yintercept = 0.9793, linetype = "dotted") +
  geom_text(aes(50,3.2724,label = "p=0.001", vjust = -0.5, hjust = 1)) +
  geom_text(aes(50,0.9793,label = "p=0.05", vjust = -0.5, hjust = 1)) +
  geom_text(aes(label = PC, vjust = -1), check_overlap = T) +
  labs(y = "Tracy Widom statistic", caption = "Tracy Widom test")


# PC loadings
p5 <- VizDimLoadings(merged.seurat, dims = c(1), reduction = "harmony") + labs(caption = "loadings for first PC")
p6 <- VizDimLoadings(merged.seurat, dims = c(2), reduction = "harmony") + labs(caption = "loadings for 2nd PC")
p7 <- VizDimLoadings(merged.seurat, dims = c(3), reduction = "harmony") + labs(caption = "loadings for 3rd PC")
p8 <- VizDimLoadings(merged.seurat, dims = c(4), reduction = "harmony") + labs(caption = "loadings for 4th PC")




ggsave(file.path(results_dir, paste0(out_stem,"_p1_master_seurat_PCA.pdf")), plot_grid(p2,p3,p5,p6,p7,p8,ncol = 3), width = 10, height = 16)

#########################################################
################### export results ######################
#########################################################




pca.meta <- merged.seurat@reductions$harmony@cell.embeddings %>% 
    		as.data.frame() %>%
    		tibble::rownames_to_column("cell.id") %>%
    		left_join(merged.seurat@meta.data %>% tibble::rownames_to_column("cell.id"), by = "cell.id") 

saveRDS(pca.meta, file.path(results_dir, paste0(out_stem,"_master_seurat_meta_table.rds")))

saveRDS(merged.seurat, file.path(objects_dir,paste0("master_seurat_post_",out_stem,".rds")))


