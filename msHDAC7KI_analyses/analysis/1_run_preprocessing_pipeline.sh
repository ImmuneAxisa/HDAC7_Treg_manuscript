#!/bin/bash

# this script will run the whole pipeline for the scRNAseq analyses. It requires R3.5


########################################################################################################
################################ First pass analysis on all data #######################################
########################################################################################################

########
## QC ##
########
Rscript A_QC/A01_create_master_seurat.R

Rscript A_QC/A02_filter_cells.R

##############################################
## dimensionality reduction and integration ##
##############################################

Rscript B_reductions/B01_PCA_harmony.R

Rscript B_reductions/B02_plot_PC_loadings.R

########################################
## clustering, visualization, markers ##
########################################

# clustering using the harmony corrected PCs
Rscript C_clustering_typing/C01_clustering.R

# clustering using the original PCs: not used
#Rscript C_clustering_typing/C02_clustering_oriPCA.R

# plots for UMAP and clustering parameter search
Rscript C_clustering_typing/C03_clustering_plots.R

# Find cluster markers, for curated cluster calls
Rscript C_clustering_typing/C04_cluster_markers.R 

####################################################################################################
################################ Second pass on T cells only #######################################
####################################################################################################

##############################################
## dimensionality reduction and integration ##
##############################################

# select the T cells and rerun PCA and harmony 
Rscript B_reductions/B03_PCA_harmony_Tcells.R

Rscript B_reductions/B04_plot_PC_loadings_Tcells.R

########################################
## clustering, visualization, markers ##
########################################

# clustering using the harmony corrected PCs
Rscript C_clustering_typing/C05_clustering_Tcells.R

# plots for UMAP and clustering parameter search
Rscript C_clustering_typing/C06_clustering_plots_Tcells.R

# Find cluster markers, for curated cluster calls
Rscript C_clustering_typing/C07_cluster_markers_Tcells.R


#######################
## DGE poisson model ##
#######################

Rscript D_DGE/D_speedglm_poisson.R


