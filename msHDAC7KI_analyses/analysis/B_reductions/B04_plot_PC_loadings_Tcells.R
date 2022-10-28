

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
library(sp)
library(harmony)



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
################# Input/output info #####################
#########################################################


# Set wd on project folder
setwd("/home/pa326/project/HDAC7_KI_EAE_10x/")

results_dir <- "./results/reductions/"
dir.create(results_dir, showWarnings = FALSE)


#########################################################
######################## plots ##########################
#########################################################

plots <- list()

# firs iteration of harmony


merged.seurat <- readRDS("./results/seurat_objects/master_seurat_post_harmony_Tcells.rds")

meta <- merged.seurat@meta.data %>% 
            tibble::rownames_to_column("cell.id")


plots[["pca.plots"]] <- merged.seurat@reductions$pca@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
		gather("PC", "loadings", -cell.id) %>%
		mutate(PC = parse_number(PC)) %>%
		dplyr::filter(PC < 11) %>%
		left_join(meta, by = "cell.id") %>%
		ggplot(aes(paste(Xp,sample), loadings, fill =Geno)) +
		    geom_boxplot() +
		    facet_grid(PC~., scales = "free_y") +
		    ggtitle("initial PCA") +
		    theme_mrl()


plots[["emulsion.harmony"]] <- merged.seurat@reductions$harmony@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
		gather("PC", "loadings", -cell.id) %>%
		mutate(PC = parse_number(PC)) %>%
		dplyr::filter(PC < 11) %>%
		left_join(meta, by = "cell.id") %>%
		ggplot(aes(paste(Xp,sample), loadings, fill =Geno)) +
		    geom_boxplot() +
		    facet_grid(PC~., scales = "free_y") +
		    ggtitle("harmony on emulsion") +
		    theme_mrl()


		    
# save plots

ggsave("./results/reductions/PC_loadings_w_integration_Tcells.pdf", 
        marrangeGrob(plots, nrow=1, ncol=1), 
        width = 30, height = 50, units = "cm")


