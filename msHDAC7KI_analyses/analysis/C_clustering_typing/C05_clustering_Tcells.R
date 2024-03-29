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
################### get seurat object ###################
#########################################################

merged.seurat <- readRDS(file.path(objects_dir,"master_seurat_post_harmony_Tcells.rds"))

cat("\n dim of seurat meta data \n")
dim(merged.seurat@meta.data)

# for testing
#merged.seurat <- merged.seurat[,1:1000]

#########################################################
########### clustering and cluster markers ##############
#########################################################

print("running clustering")

# clustering function

run_clustering <- function(seurat_object){

    print("find neighbors")
    for(k in c(20)) {
    
            seurat_object <- FindNeighbors(seurat_object,
                                            dims = 1:30, 
                                            k.param = k,
                                            reduction = "harmony")
    
    
    print("find clusters")
        for(algo in c(1)) {
                        
                        
                        # Loop with different resolutions to define clusters
                        
                        for(res in c(0.3, 0.4, 0.5)) {
                        
                        
                            seurat_object <- FindClusters(seurat_object, 
                                                            resolution = res,
                                                            algorithm = algo)
                            clusters <- seurat_object[["seurat_clusters"]][,1]
                            seurat_object[[paste0("SNNk",k,"_algo",algo,"_res",res)]] <- clusters
                            }
                            
                        }
                        
    }


    return(seurat_object)
}


merged.seurat <- run_clustering(merged.seurat)


#############################################################
###### run subclustering on low resolution communities ######
#############################################################


# Select a cluster iteration
cluster_assignments <- merged.seurat@meta.data$SNNk20_algo1_res0.3

# run subclustering on each cluster
subclustering <- lapply(unique(cluster_assignments), function(cluster) {


    # subset to given cluster
    tmp_seurat <- merged.seurat[, cluster_assignments == cluster]
    
    # rerun SNN graph building
    tmp_seurat <- FindNeighbors(tmp_seurat,
                                    dims = 1:30, 
                                    k.param = 20,
                                    reduction = "harmony")

    # find clusters at different resolutions
    for(res in c(0.05,0.07,0.1,0.15,0.2)) {
                        
        cat("\n\n\n subclustering resolution:",res,"\n\n\n")                
        tmp_seurat <- FindClusters(tmp_seurat, 
                                        resolution = res,
                                        algorithm = 1)
        
        # extract clustering results                                
        subclusters <- tmp_seurat[["seurat_clusters"]][,1]
        
        # rename it and add back to metadata
        tmp_seurat[[paste0("Sub_res",res)]] <- paste0(cluster,"_",subclusters)
    }
    
    
    # only keep the columns with the results of the subclustering, and return the df on subcluster calls
    tmp_meta <- tmp_seurat@meta.data %>%
        tibble::rownames_to_column("cell.id") %>%
        dplyr::select(matches("Sub_res|cell.id"))
    
})

# merge all the dataframe corresponding to each cluster
subclustering <- do.call(rbind, subclustering)


# add the results back to the main object
meta <- merged.seurat@meta.data %>% 
    tibble::rownames_to_column("cell.id") %>%
    left_join(subclustering) %>%
    tibble::column_to_rownames("cell.id")


merged.seurat@meta.data <- meta



#########################################################
############## tSNE and UMAP reductions #################
#########################################################

print("tSNE")
merged.seurat <- RunTSNE(merged.seurat, reduction = "harmony", dims = 1:30)

reductions.meta <- merged.seurat@reductions$tsne@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
    	left_join(merged.seurat@meta.data %>% tibble::rownames_to_column("cell.id"), by = "cell.id") 


print("running UMAP iterations")

for(n in c(50)) {
    for(min.dist in c(0.01,0.3,0.5)) {
        for(spread in c(5,10)) {
            try(if(T){ 
                it <- paste0("UMAPn",n,"xMd",min.dist,"xS", spread)
                message(paste0("running ", it))
                merged.seurat <- RunUMAP(merged.seurat, 
                                        dims = 1:30, 
                                        reduction = "harmony", 
                                        n.neighbors = n,
                                        min.dist = min.dist, 
                                        spread = spread,
                                        reduction.key = it)
                                           
                reductions.meta <- merged.seurat@reductions$umap@cell.embeddings %>% 
                		as.data.frame() %>%
                		tibble::rownames_to_column("cell.id") %>%
                		left_join(reductions.meta, by = "cell.id") 
            })    		
        }
    }                        
}





#########################################################
#################### Save outputs #######################
#########################################################

saveRDS(reductions.meta, file.path(results_dir,"master_seurat_Tcells_low_res_clustering_table_harmony.rds"))


saveRDS(merged.seurat, file.path(objects_dir,"master_seurat_harmony_Tcells_low_res_clustered.rds"))
