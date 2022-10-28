#################################################
################### Set up ######################
#################################################



setwd("/gpfs/ycga/home/pa326/project/HDAC7_KI_EAE_10x")

results_dir <- "results/master"
dir.create(results_dir, showWarnings = FALSE)
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




#################################################
################ Import data ####################
#################################################

# Set path to cell ranger outputs  

path <- "/gpfs/ycga/scratch60/hafler/pa326/HDAC7_KI_EAE_10x/cellranger/"


samples <- list.files(path)

## Import expression matrix as a list of seurat objects
seurat_raw <- lapply(samples, function(x) {
                h5.path <- paste0(path,x,"/filtered_feature_bc_matrix.h5")
                Seurat_object <- CreateSeuratObject(Read10X_h5(h5.path))
                return(Seurat_object)
})

# name list elements with sample names
names(seurat_raw) <- gsub("_MMT_cellranger","", samples)



#################################################
############### merge objects ###################
#################################################

# merge all seurat objects and add a suffix to each cell barcode with the sample names it belongs to
merged.seurat <- merge(seurat_raw[[1]], seurat_raw[-1], add.cell.ids = names(seurat_raw))


#################################################
################ add metadata ###################
#################################################

# Add gene metadata

meta.features <- merged.seurat@assays[["RNA"]]@meta.features %>%
    tibble::rownames_to_column("gene_name") %>%
    left_join(readRDS("./data/gene_annotations_positions.rds")) %>%
    tibble::column_to_rownames("gene_name")

merged.seurat@assays[["RNA"]]@meta.features <- meta.features


# add sample level metadata

sample_meta <- read_csv("data/sample_meta.csv")


raw_meta <- merged.seurat@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  separate(cell.id, c("sample", "barcode"), sep = "_", remove = F) %>%
  left_join(sample_meta) %>%
  dplyr::select(-barcode) 

#################################################
########## filter out unwanted samples ##########
#################################################

# remove samples using metadata
cat("\n n cells before filtering:", nrow(raw_meta))

raw_meta <-  raw_meta
    
cat("\n n cells after filtering:", nrow(raw_meta))

# subset object with remaining cell.ids in meta
merged.seurat <- merged.seurat[, pull(raw_meta, "cell.id")]


cat("\n n cells in seurat object after filtering", ncol(merged.seurat))



# add updated meta back into the object
merged.seurat@meta.data <- raw_meta %>%
  tibble::column_to_rownames("cell.id")




####################################################
################ calc QC metrics ###################
####################################################

merged.seurat[["percent.mt"]] <- PercentageFeatureSet(merged.seurat, pattern = "^mt-")
merged.seurat[["percent.ribo"]] <- PercentageFeatureSet(merged.seurat, pattern = "^Rp[rl]")


meta.features <- tibble::rownames_to_column(meta.features,"gene_name")

merged.seurat[["percent.TCRab"]] <- PercentageFeatureSet(merged.seurat, features = dplyr::filter(meta.features, gene_category == "TCRab") %>% pull(gene_name))
merged.seurat[["percent.TCRgd"]] <- PercentageFeatureSet(merged.seurat, features = dplyr::filter(meta.features, gene_category == "TCRgd") %>% pull(gene_name))
merged.seurat[["percent.BCR"]] <- PercentageFeatureSet(merged.seurat, features = dplyr::filter(meta.features, gene_category == "BCR") %>% pull(gene_name))




#################################################
################ save to file ###################
#################################################


saveRDS(merged.seurat, file.path(results_dir,"seurat_objects","master_seurat.rds"))

saveRDS(merged.seurat@meta.data, file.path(results_dir,"meta","master_meta.rds"))


