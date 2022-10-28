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
#library(hdf5r)
library(ggsci)
#library(harmony)
library(MASS)
library(reshape2)
library(forcats)



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
#################### get seurat meta ####################
#########################################################

meta_table <- readRDS(file.path(results_dir,"master_seurat_Tcells_low_res_clustering_table_harmony.rds"))


#########################################################
########### plotting defaults and functions #############
#########################################################



options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

palette <- c(pal_rickandmorty()(12),pal_jco()(10), pal_futurama()(11), pal_tron()(7))

##-----------------------------------
## 2D density by probability function


density2D <- function(mv, dim1, dim2, groups = NULL, smoothing = 1) {
  
  smoothing = 1.2
  
  # core function
  .density2D <- function(mv, dim1, dim2) {
    #print(mv[,dim1])
    mv.kde <- kde2d(mv[,dim1], 
                    mv[,dim2], 
                    h = c(bandwidth.nrd(mv[,dim1])*smoothing,
                          bandwidth.nrd(mv[,dim2])*smoothing),
                    n = 100)
    
    
    dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
    dy <- diff(mv.kde$y[1:2])
    sz <- sort(mv.kde$z)
    c1 <- cumsum(sz) * dx * dy
    
    dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
    dc <- melt(mv.kde$z)
    dc$prob <- approx(sz,1-c1,dc$value)$y
    return(dc)
  }
  
  
  # split by group if specified
  if(is.null(groups)) {
    dc <- .density2D(mv, dim1, dim2)
  } else {
    #print(unique(pull(mv,groups)))
    
    # Iterate over each group
    dc.list <- lapply(unique(pull(mv,groups)), function(group) {
      #print(group)
      # Subset the df for each group
      mv.subset <- mv %>%
        dplyr::filter(!!rlang::sym(groups) == group)
      #print(mv.subset)
      
      # Run density2D on the subset
      mv.subset <- .density2D(mv.subset, dim1, dim2) %>%
        mutate(!!quo_name(groups) := group)
    })
    dc <- do.call(rbind, dc.list)
  }
  
  # Rename columns to original name
  
  dc <- dc %>%
    dplyr::rename(!!quo_name(dim1)  := Var1, !!quo_name(dim2)  := Var2)

}


##-----------------------------------
## 






###########################################
########### rearrange results #############
###########################################


# gather UMAP low dimensionality embeddings coordinates

UMAP <- meta_table %>%
  dplyr::select(matches("UMAP|cell.id")) %>%
  gather("embed", "coord", -cell.id) %>%
  separate(embed, c("iteration", "dim"), sep = "_") %>%
  mutate(dim = paste0("dim_",dim)) %>%
  spread(dim,coord) %>%
  separate(iteration, c("knn", "min.dist","spread"), sep = "x") %>%
  mutate(min.dist = gsub("Md0", "Md0.", min.dist)) %>%
  mutate(min.dist = fct_reorder(min.dist, 
                                as.numeric(gsub("Md0", "0", min.dist))
                                )) %>%
  mutate(spread = fct_reorder(spread, 
                                as.numeric(gsub("S", "", spread))
                                ))

# gather main cluster calls

clusters <- meta_table %>%
  dplyr::select(matches("SNN|cell.id")) %>%
  gather("cluster_param", "cluster_call", -cell.id) %>%
  separate(cluster_param, c("SNNk", "algo","res"), sep = "_")



# gather sub cluster calls

subclusters <- meta_table %>%
  dplyr::select(matches("Sub_res|cell.id")) %>%
  gather("Sub_res", "subcluster_call", -cell.id)
  
###########################################
############# which plots #################
###########################################

UMAP_topo <- T
clusters_topo <- T
clusters_concentration <- T
subclusters_topo <- T

results <- list()  
  


###########################################
############# UMAP results ################
###########################################

if(UMAP_topo) {

UMAP_topo <- lapply(unique(UMAP$knn),
                    function(knn_index) {
                        
                        plot <- UMAP %>%
                            dplyr::filter(knn == knn_index) %>%
                            #mutate(dim_2 = ifelse(knn == "UMAPn50", -dim_2, dim_2)) %>% # invert dim 2 for knn = 50 as it seems flipped
                            ggplot(aes(dim_1,
                                       dim_2)) +
                            theme_mrl() +
                            geom_hex(aes(color = ..count..), bins = 200) +
                            scale_color_viridis(trans = "log2") +
                            scale_fill_viridis(trans = "log2") +
                            facet_wrap(min.dist~spread, scales = "free") +
                            theme(axis.line.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.line.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank()) +
                            ggtitle(paste("knn =",knn_index))
 
                        
                    })


results[["UMAP_topo"]] <- c(list(ggplot() + 
                                    ggtitle("\n\n\n\n UMAP topology") + 
                                    theme_void() +
                                    theme(plot.title =  element_text(size = 30, hjust=0.5, vjust=2, face="bold"))),
                                UMAP_topo)

}

##################################
########## clustering ############
##################################


combinations <- expand.grid(knn = unique(UMAP$knn),
                            cluster = grep("SNNk", colnames(meta_table), value = T),
                            stringsAsFactors = F) %>%
                        purrr::transpose()


if(clusters_topo) {

clusters_topo <- lapply(combinations,
                    function(combinations_index) {
                    
                        clusters <- meta_table %>%
                            dplyr::select("cell.id", combinations_index[["cluster"]]) %>%
                            dplyr::rename("clustering" = combinations_index[["cluster"]])
                        
                        plot <- UMAP %>%
                            dplyr::filter(knn == combinations_index[["knn"]]) %>%
                            left_join(clusters) %>%
                            #mutate(dim_2 = ifelse(knn == "UMAPn50", dim_2, -dim_2)) %>% # invert dim 2 for knn = 50 as it seems flipped
                            ggplot(aes(dim_1,
                                       dim_2,
                                       color = clustering,
                                       fill = clustering)) +
                            stat_density2d(bins =5, 
                                            contour_var = "ndensity", 
                                            alpha = 0.3,
                                            geom = "polygon") +
                            facet_wrap(min.dist~spread, scales = "free") +
                            theme_mrl() +
                            theme(axis.line.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.line.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank()) +
                            ggtitle(paste(combinations_index[["knn"]],
                                            "\n",
                                            combinations_index[["cluster"]]))
                        
                    })
                    

results[["clusters_topo"]] <- c(list(ggplot() + 
                                        ggtitle("\n\n\n\n clusters_topo") + 
                                        theme_void() +
                                        theme(plot.title =  element_text(size = 30, hjust=0.5, vjust=2, face="bold"))),
                                clusters_topo)

}





if(clusters_concentration) {

                    
                    
clusters_concentration <- lapply(combinations,
                    function(combinations_index) {
                    
                        clusters <- meta_table %>%
                            dplyr::select("cell.id", combinations_index[["cluster"]]) %>%
                            dplyr::rename("clustering" = combinations_index[["cluster"]])
                        
                        tmp <- UMAP %>%
                            dplyr::filter(knn == combinations_index[["knn"]] & grepl("5",spread) & min.dist ==  "Md0.01") 
                            
                        
                            #mutate(dim_2 = ifelse(knn == "UMAPn50", dim_2, -dim_2)) %>% # invert dim 2 for knn = 50 as it seems flipped
                        plot <- tmp %>%
                            left_join(clusters) %>%
                            ggplot(aes(dim_1,
                                       dim_2)) +
                            stat_bin_hex(data = tmp, 
                                         color = "grey",
                                         fill = "grey",
                                         bins = 50) +
                            stat_bin_hex(aes(alpha = ..ndensity..),
                                         fill = "orange",
                                         color = "black",
                                         lwd = 0.1,
                                         bins = 50) +
                            facet_wrap(.~clustering) +
                            theme_mrl() +
                            theme(axis.line.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.line.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank()) +
                            ggtitle(paste(combinations_index[["knn"]],
                                            "\n",
                                            combinations_index[["cluster"]]))
                        
                    })

results[["clusters_concentration"]] <- c(list(ggplot() + 
                                                ggtitle("\n\n\n\n clusters_concentration") + 
                                                theme_void() +
                                                theme(plot.title =  element_text(size = 30, hjust=0.5, vjust=2, face="bold"))),
                                        clusters_concentration)


}
                            
##################################
########## subclustering  ########
##################################


combinations <- expand.grid(knn = unique(UMAP$knn),
                            cluster = grep("Sub_res", colnames(meta_table), value = T),
                            stringsAsFactors = F) %>%
                        purrr::transpose()
                        
if(subclusters_topo) {

subclusters_topo <- lapply(combinations,
                    function(combinations_index) {
                    
                        clusters <- meta_table %>%
                            dplyr::select("cell.id", combinations_index[["cluster"]]) %>%
                            dplyr::rename("clustering" = combinations_index[["cluster"]])
                        
                        plot <- UMAP %>%
                            dplyr::filter(knn == combinations_index[["knn"]]) %>%
                            left_join(clusters) %>%
                            #mutate(dim_2 = ifelse(knn == "UMAPn50", dim_2, -dim_2)) %>% # invert dim 2 for knn = 50 as it seems flipped
                            ggplot(aes(dim_1,
                                       dim_2,
                                       color = clustering,
                                       fill = clustering)) +
                            stat_density2d(bins =5, 
                                            contour_var = "ndensity", 
                                            alpha = 0.3,
                                            geom = "polygon") +
                            facet_wrap(min.dist~spread, scales = "free") +
                            theme_mrl() +
                            theme(axis.line.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.line.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank()) +
                            ggtitle(paste(combinations_index[["knn"]],
                                            "\n",
                                            combinations_index[["cluster"]]))
                        
                    })
                    

results[["subclusters_topo"]] <- c(list(ggplot() + 
                                        ggtitle("\n\n\n\n subclusters_topo") + 
                                        theme_void() +
                                        theme(plot.title =  element_text(size = 30, hjust=0.5, vjust=2, face="bold"))),
                                subclusters_topo)
                                
                                
}

##############################
########## export  ###########
##############################


results <- do.call(c,results)



ggsave(file.path(results_dir, "low_res_clustering_plots_Tcells_harmony.pdf"), 
        marrangeGrob(results, nrow=1, ncol=1), 
        width = 40, height = 30, units = "cm")





