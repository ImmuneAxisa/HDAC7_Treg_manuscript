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
library(rhdf5)


#################################################
################ Import data ####################
#################################################

# read file
merged.seurat <- readRDS(file.path(results_dir,"seurat_objects","master_seurat.rds"))

# get counts
count_matrix <- as.matrix(merged.seurat@assays$RNA@counts)

rm(merged.seurat)
gc()

# save h5
h5createFile(file.path(results_dir,"raw_counts.h5"))

h5createDataset(file.path(results_dir,"raw_counts.h5"), 'raw_counts', dim(count_matrix), chunk = c(1000, 1000))
h5write(count_matrix, file.path(results_dir,"raw_counts.h5"), "raw_counts")

#h5writeDataset(count_matrix, file.path(results_dir,"raw_counts.h5"), "raw_counts")
