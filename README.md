# Introduction

This repository contains the code associated with the manuscript "A Multiple Sclerosis Protective Coding Variant Reveals an Essential Role for HDAC7 in Regulatory T cells". This README file describes the organisation of the scripts related to the analyses presented in the manuscript. This repository has also been deposited on figshare.com (doi:10.6084/m9.figshare.21539098), along with various intermediate results file that might be useful to users who wish to reuse the data and analyses in this project. The github repository only contains the scripts, but the figshare.com repository contains the intermediate results file as well. For ease of use, the full directory and file tree structure deposited on figshare.com is presented at the end of this document.

# Disclaimer
* Throughout this project, We use the `SLURM` job scheduling system to run some of the "heavy" analyses on the Yale University Ruddle cluster.
* The code and additional resources on figshare.com are provided as-is. We have validated that this runs without errors on the Ruddle cluster and on local MacOS systems, but do not guarantee that this will be the case in different computing environments.
* The "light" analysis portion of this code is written with relative paths that have root project directory structure within this repository (detailed below) but the paths in the scripts run on Ruddle will need to be adjusted to meet you user names and path structures.
* All analyses run in R were run with `R3.5` (unless stated otherwise in the script header). We do not guarantee the scripts will run with more recent versions of R

# Content

The repository is split into 2 main parts, each of which serve as project root for the scripts contained in them.
* human_datasets_analyses: contains all human transcriptmic data analysed, generated in this manuscript or reanalyses of publicly available datasets.
* msHDAC7KI_analyses: contains all analyses related to the in vivo EAE model using HDAC7 R150KI mice.

## human_datasets_analyses

This directory is the project root for all scripts in the subdirectories.

### src

* `bulk.rna.align.sh`: alignment and quantitations from the original Fastq files. This script is designed to be run on the Ruddle cluster with the `SLURM` scheduler distributing one job per sample. To do so, the scripts iterates through a sample table specifying the name of the fastq files associated with each sample. The sample table is expected to contain at least the following columns: Sample_Number, Sample_Name, RNA_FastQ1, RNA_FastQ2.
* `bulk_wrappers.R` and `gseaplot2_PPAmod.R` contain custom functions for differential gene expression analyses in R. They are sourced by the script in the analysis folder.


### huHDAC7_coexpression

Co-expression analyses on human Treg cells, highlighting the main module containing HDAC7, related to Figure 1. Using this module, we run enrichment analyses against public databases using enrichr, and against a previously published single cell RNAseq dataset of human Treg cells in skin and blood.

### huHDAC7_KD

Analyses related to the knock-down of HDAC7 in human Treg cells presented in Figure 2 and supplementary figure 2. We run differential expression analyses with DESeq2 and then run GSEA on custom FOXO geneset generated by reanalysing publicly available mouse Treg datasets KO for FOXO1/3.


### huHDAC7_variant

Analyses related to the knock-down of HDAC7 in human Treg cells presented in supplementary figure 6. We run differential expression analyses with DESeq2 and then run GSEA on hallmarks genesets, and geneset from the HDAC7 R150H KI mouse Treg infiltrating the brain during EAE (presented below).

## msHDAC7KI_Analyses

Single cell RNAseq of brain infiltrating T cells at the peak of EAE in HDAC7 R150H KI mice vs WT controls, presented in Figure 5, supplementary figure 10, and supplementary figure 11. This directory is the project root for all analyses in the subdirectories.

### src

Contains the script to run cellranger on the Fastq files and generated a count matrix for all droplets. This is a simple adaptation of the manual instructions available on the 10x genomics website, to be run on the Ruddle cluster.

### analysis: 

It is split into 2 portions:

* heavy analyses run on the Ruddle cluster are located in the subfolders ranging from A to D, and can be called in order by `1_run_preprocessing_pipeline.sh`. The scripts called assume the directory structure on Ruddle and hence paths need to be adjusted for these to run properly. The pipeline is structured as follows:

  * gather results from the cell ranger output and create a seurat object (A_01)
  * remove bad quality cells and potential doublets (A_02)
  * Run PCA and integrate across emulsions using `harmony` to generate corrected PC loadings (B_01) and plot PC loadings (B_02)
  * run Louvain clustering and UMAP with a parameter search (C_01) and plot the results (C03)
  * Based on plots curate clusters and run the `presto` find markers function to define cluster markers (C04)
  * Select T cell clusters and rerun PCA and `harmony` (B_03) and plot PC loadings (B_04).
  * Rerun Louvain clustering and UMAP with a parameter search (C_05) and plot the results (C06)
  * Based on plots curate clusters and run the `presto` find markers function to define cluster markers (C07)
  * Run DGE using a poisson model with the `speedglm` implementation (D). This step is heavy to run and we use the `rslurm` package to parallelize the computation.
  
* `2_MASTER_clustering_typing.R` gathers the results from the preprocessing pipeline and plots the results of dimensionality reductions, clustering, cell typing, marker genes expression.

* `3_MASTER_DGE.R` gather the results from the poisson model run in D. And plot the results.



# Full directory tree containing scripts and intermediate results files, tables and figures

```bash
├── README.md
├── human_datasets_analyses
│   ├── huHDAC7_KD
│   │   ├── analysis
│   │   │   ├── DGE.R                         # DESeq2 differential gene expression analysis
│   │   │   └── GSEA.R                        # GSEA on DESeq2 outputs, using the FOXO KO genesets
│   │   ├── data
│   │   │   ├── GEO_table_STM_huKD.csv        # sample table
│   │   │   ├── foxo1.3.TregKO.txt            # GEO2R results table for FOXO1/3 DKO mouse Treg
│   │   │   ├── foxo1.TregKO.txt              # GEO2R results table for FOXO1 KO Treg
│   │   │   └── lost-retrieved-from-affy.tsv  # probe id to gene correspondance
│   │   └── results
│   │       ├── DGE
│   │       │   └── huHDAC7_KD_STM.csv        # DESeq2 results table
│   │       ├── figs
│   │       │   ├── Foxos_genesets.pdf        
│   │       │   ├── HDAC7exp.pdf
│   │       │   ├── Heatmap_STM.pdf
│   │       │   ├── MAplot.pdf
│   │       │   ├── batchedCorrectedPCA.pdf
│   │       │   ├── gseaplots.pdf
│   │       │   └── gseaplots2.pdf
│   │       ├── raw_counts.h5                 # raw count matrix uploaded to GEO
│   │       └── rna                           # outputs from running bulk.rna.align.sh on the fastq
│   │           └── quantitation
│   │               └── rsem
│   │                   ├── TS12.genes.results
│   │                   ├── TS13.genes.results
│   │                   ├── TS2.genes.results
│   │                   ├── TS3.genes.results
│   │                   ├── TS7.genes.results
│   │                   └── TS8.genes.results
│   ├── huHDAC7_coexpression
│   │   ├── analysis
│   │   │   ├── TeichmannTreg_reanalysis.R                  # scRNAseq of huTreg reanalysis, scoring HDAC7 module
│   │   │   └── wgcna_and_enrichr.R                         # coexpression networks in Tregs and enrichment on HDAC7 module
│   │   └── results
│   │       ├── HDAC7moduleTreg.txt                         # module 1 gene name list containing HDAC7
│   │       ├── HDAC7moduleTreg_wID.txt                     # module 1 gene table containing HDAC7 with gene name and IDs
│   │       ├── TeichTreg
│   │       │   └── Teichman.Treg.skinblood.seurat.rds      # exported seurat object after reanalysis
│   │       ├── enrichr                                     # results from enrichr online portal analysis
│   │       │   ├── GO_Cellular_Component_2018_table.txt 
│   │       │   └── Transcription_Factor_PPIs_table.txt
│   │       ├── figs
│   │       │   ├── GO_CC.pdf
│   │       │   ├── HDAC7moduleScore_STM.pdf
│   │       │   ├── HDAC7moduleScore_cluster_STM.pdf
│   │       │   ├── Module_colors.pdf
│   │       │   ├── PPIenrichment.pdf
│   │       │   ├── WGCNAdendrogramTregsLarge.pdf
│   │       │   ├── clusterUMAP_STM.pdf
│   │       │   ├── cluster_tissue_plot.pdf
│   │       │   ├── dotplot.pdf
│   │       │   └── tissueUMAP_STM.pdf
│   │       └── treg.rsem.filtered.counts.txt
│   ├── huHDAC7_variant
│   │   ├── analysis
│   │   │   ├── DGE.R                                                         # DESeq2 differential gene expression analysis
│   │   │   └── GSEA.R                                                        # GSEA on hallmarks and msHDAC7 R150H KI Treg set
│   │   ├── data
│   │   │   ├── GEO_table_STM_variant.csv                                     # sample table
│   │   │   ├── Homo_sapiens.GRCh38.88.primary_assembly.gene_symbols.txt      # gene annotations
│   │   │   └── h.all.v6.2.symbols.gmt                                        # hallmark genesets
│   │   └── results
│   │       ├── DGE                                                           # DESeq2 results tables
│   │       │   ├── degDESeqMut.txt
│   │       │   ├── degDESeqVariant.txt
│   │       │   └── degDESeqWT.txt
│   │       ├── figs
│   │       │   ├── HDAC7exp.pdf
│   │       │   ├── batchedCorrectedPCA.pdf
│   │       │   ├── diagonal_plot_2.pdf
│   │       │   ├── gseaplots.pdf
│   │       │   ├── gseaplots2.pdf
│   │       │   ├── heatmap.pdf
│   │       │   └── huVar_vs_msVar.pdf
│   │       ├── raw_counts.h5                                                 # raw count matrix uploaded to GEO
│   │       └── rna                                                           # outputs from running bulk.rna.align.sh on the fastq
│   │           └── quantitation
│   │               └── rsem
│   │                   ├── P13.genes.results
│   │                   ├── P14.genes.results
│   │                   ├── P15.genes.results
│   │                   ├── P22.genes.results
│   │                   ├── P23.genes.results
│   │                   ├── P24.genes.results
│   │                   ├── P5.genes.results
│   │                   ├── P6.genes.results
│   │                   ├── PP1.genes.results
│   │                   ├── PP2R.genes.results
│   │                   ├── PP3R.genes.results
│   │                   ├── PP4R.genes.results
│   │                   ├── PP5R.genes.results
│   │                   ├── PP6R.genes.results
│   │                   ├── PP7.genes.results
│   │                   ├── PP8.genes.results
│   │                   └── PP9.genes.results
│   └── src             
│       ├── bulk.rna.align.sh     # alignment and quantitation script. Run on Yale Ruddle cluster with SLURM scheduler
│       ├── bulk_wrappers.R       # custom function wrappers for DGE and GSEA
│       └── gseaplot2_PPAmod.R    # custom GSEA plot function
└── msHDAC7KI_analyses
    ├── analysis
    │   ├── 1_run_preprocessing_pipeline.sh               # preprocessing pipeline. Run on Yale Ruddle cluster with SLURM scheduler 
    │   ├── 2_MASTER_clusturing_typing.R                  # Plots results for dim reduction, clustering, cluster markers
    │   ├── 3_MASTER_DGE.R                                # runs poisson model DGE on clusters of interest: Th17 and Treg clusters
    │   ├── A_QC                                          # Pipeline first step: create object and QC emulsions
    │   │   ├── A01_create_master_seurat.R
    │   │   └── A02_filter_cells.R
    │   ├── B_reductions                                  # Pipeline second step: PCA, harmony.
    │   │   ├── B01_PCA_harmony.R
    │   │   ├── B02_plot_PC_loadings.R
    │   │   ├── B03_PCA_harmony_Tcells.R
    │   │   └── B04_plot_PC_loadings_Tcells.R
    │   ├── C_clustering_typing                           # Pipeline 3rd step: clustering, UMAPs, cluster markers
    │   │   ├── C01_clustering.R
    │   │   ├── C02_clustering_oriPCA.R
    │   │   ├── C03_clustering_plots.R
    │   │   ├── C04_cluster_markers.R
    │   │   ├── C05_clustering_Tcells.R
    │   │   ├── C06_clustering_plots_Tcells.R
    │   │   └── C07_cluster_markers_Tcells.R
    │   ├── D_DGE                                         # Pipeline 4th step: differential gene expression
    │   │   └── D_speedglm_poisson.R
    │   └── GEO_matrix_generation.R                       # Generate raw matrix for GEO upload
    ├── data
    │   ├── gene_annotations_positions.rds                # gene id to gene name correspondance. With link to Seurat names (eg .1 suffixes)
    │   └── sample_meta.csv                               # sample table
    ├── figures
    │   └── master
    │       ├── DGE                                       # figures for DGE results on each clusters. pdf for paper figures.
    │       │   ├── all 0_0ma_plot_ashr_fdr.png
    │       │   ├── all 0_0p_val_histogram.png
    │       │   ├── all 0_0volcano_plot.png
    │       │   ├── all 0_0volcano_plot_ashr_fdr.png
    │       │   ├── all 0_1ma_plot_ashr_fdr.png
    │       │   ├── all 0_1p_val_histogram.png
    │       │   ├── all 0_1volcano_plot.png
    │       │   ├── all 0_1volcano_plot_ashr_fdr.png
    │       │   ├── all 4_0ma_plot_ashr_fdr.png
    │       │   ├── all 4_0p_val_histogram.png
    │       │   ├── all 4_0volcano_plot.png
    │       │   ├── all 4_0volcano_plot_ashr_fdr.png
    │       │   ├── all 6_0ma_plot_ashr_fdr.png
    │       │   ├── all 6_0p_val_histogram.png
    │       │   ├── all 6_0volcano_plot.png
    │       │   ├── all 6_0volcano_plot_ashr_fdr.png
    │       │   ├── all 6_1ma_plot_ashr_fdr.png
    │       │   ├── all 6_1p_val_histogram.png
    │       │   ├── all 6_1volcano_plot.png
    │       │   ├── all 6_1volcano_plot_ashr_fdr.png
    │       │   ├── paper_MAplot_Teff.pdf
    │       │   ├── paper_MAplot_Treg.pdf
    │       │   ├── paper_MAplot_Tresting.pdf
    │       │   └── paper_MAplot_central.pdf
    │       ├── clusters                                  # figure for clustering and UMAP results. main cell types and fine res cell types
    │       │   ├── fine.png
    │       │   ├── fine_legend.pdf
    │       │   ├── freq.pdf
    │       │   ├── main.png
    │       │   └── main_legend.pdf
    │       └── marker_genes                              # marker gene expression overlayed on UMAPs
    │           ├── gene_Ass1_mean.png
    │           ├── gene_Bcl2_mean.png
    │           ├── gene_Bcl2a1a_mean.png
    │           ├── gene_Bcl2a1b_mean.png
    │           ├── gene_Bcl2a1d_mean.png
    │           ├── gene_Bcl2l11_mean.png
    │           ├── gene_Car2_mean.png
    │           ├── gene_Ccl1_mean.png
    │           ├── gene_Ccl4_mean.png
    │           ├── gene_Ccl5_mean.png
    │           ├── gene_Ccr2_mean.png
    │           ├── gene_Ccr7_mean.png
    │           ├── gene_Ccr8_mean.png
    │           ├── gene_Cd164_mean.png
    │           ├── gene_Cd226_mean.png
    │           ├── gene_Cd244a_mean.png
    │           ├── gene_Cd247_mean.png
    │           ├── gene_Cd274_mean.png
    │           ├── gene_Cd27_mean.png
    │           ├── gene_Cd28_mean.png
    │           ├── gene_Cd2_mean.png
    │           ├── gene_Cd3e_mean.png
    │           ├── gene_Cd40lg_mean.png
    │           ├── gene_Cd44_mean.png
    │           ├── gene_Cd47_mean.png
    │           ├── gene_Cd48_mean.png
    │           ├── gene_Cd4_mean.png
    │           ├── gene_Cd52_mean.png
    │           ├── gene_Cd53_mean.png
    │           ├── gene_Cd5_mean.png
    │           ├── gene_Cd69_mean.png
    │           ├── gene_Cd6_mean.png
    │           ├── gene_Cd72_mean.png
    │           ├── gene_Cd74_mean.png
    │           ├── gene_Cd7_mean.png
    │           ├── gene_Cd8a_mean.png
    │           ├── gene_Cd8b1_mean.png
    │           ├── gene_Cd9_mean.png
    │           ├── gene_Csf2_mean.png
    │           ├── gene_Ctla2a_mean.png
    │           ├── gene_Ctla4_mean.png
    │           ├── gene_Cxcr3_mean.png
    │           ├── gene_Cxcr6_mean.png
    │           ├── gene_Eomes_mean.png
    │           ├── gene_Fasl_mean.png
    │           ├── gene_Fcer1g_mean.png
    │           ├── gene_Foxp3_mean.png
    │           ├── gene_Gzma_mean.png
    │           ├── gene_Gzmb_mean.png
    │           ├── gene_Gzmk_mean.png
    │           ├── gene_Hif1a_mean.png
    │           ├── gene_Icos_mean.png
    │           ├── gene_Ifng_mean.png
    │           ├── gene_Ifngr1_mean.png
    │           ├── gene_Ifngr2_mean.png
    │           ├── gene_Ikzf2_mean.png
    │           ├── gene_Il10ra_mean.png
    │           ├── gene_Il12rb2_mean.png
    │           ├── gene_Il17a_mean.png
    │           ├── gene_Il1r2_mean.png
    │           ├── gene_Il2ra_mean.png
    │           ├── gene_Il2rb_mean.png
    │           ├── gene_Il2rg_mean.png
    │           ├── gene_Il7r_mean.png
    │           ├── gene_Itgb7_mean.png
    │           ├── gene_Junb_mean.png
    │           ├── gene_Jund_mean.png
    │           ├── gene_Klf2_mean.png
    │           ├── gene_Klrb1c_mean.png
    │           ├── gene_Klrc1_mean.png
    │           ├── gene_Klrd1_mean.png
    │           ├── gene_Lag3_mean.png
    │           ├── gene_Lef1_mean.png
    │           ├── gene_Lgals1_mean.png
    │           ├── gene_Lmo4_mean.png
    │           ├── gene_Ly6a_mean.png
    │           ├── gene_Ly6c2_mean.png
    │           ├── gene_Ncr1_mean.png
    │           ├── gene_Nfkbia_mean.png
    │           ├── gene_Nkg7_mean.png
    │           ├── gene_Nr4a1_mean.png
    │           ├── gene_Pdcd1_mean.png
    │           ├── gene_Pdcd1lg2_mean.png
    │           ├── gene_Ppp1r3b_mean.png
    │           ├── gene_Rorc_mean.png
    │           ├── gene_S1pr1_mean.png
    │           ├── gene_Sell_mean.png
    │           ├── gene_Tbx21_mean.png
    │           ├── gene_Tcf7_mean.png
    │           ├── gene_Tgfb1_mean.png
    │           ├── gene_Tigit_mean.png
    │           ├── gene_Tmem176b_mean.png
    │           ├── gene_Trdv2-2_mean.png
    │           ├── gene_Trdv4_mean.png
    │           ├── gene_Xcl1_mean.png
    │           └── gene_Zbtb16_mean.png
    ├── results                 
    │   ├── cellranger                                      # results from cellranger count outputs. Only h5 files and summary stats kept
    │   │   ├── KI6_MMT_cellranger
    │   │   │   ├── filtered_feature_bc_matrix.h5
    │   │   │   ├── metrics_summary.csv
    │   │   │   ├── molecule_info.h5
    │   │   │   ├── raw_feature_bc_matrix.h5
    │   │   │   └── web_summary.html
    │   │   ├── KI7_MMT_cellranger
    │   │   │   ├── filtered_feature_bc_matrix.h5
    │   │   │   ├── metrics_summary.csv
    │   │   │   ├── molecule_info.h5
    │   │   │   ├── raw_feature_bc_matrix.h5
    │   │   │   └── web_summary.html
    │   │   ├── KI8_MMT_cellranger
    │   │   │   ├── filtered_feature_bc_matrix.h5
    │   │   │   ├── metrics_summary.csv
    │   │   │   ├── molecule_info.h5
    │   │   │   ├── raw_feature_bc_matrix.h5
    │   │   │   └── web_summary.html
    │   │   ├── KI_MMT_cellranger
    │   │   │   ├── filtered_feature_bc_matrix.h5
    │   │   │   ├── metrics_summary.csv
    │   │   │   ├── molecule_info.h5
    │   │   │   ├── raw_feature_bc_matrix.h5
    │   │   │   └── web_summary.html
    │   │   ├── WT6_MMT_cellranger
    │   │   │   ├── filtered_feature_bc_matrix.h5
    │   │   │   ├── metrics_summary.csv
    │   │   │   ├── molecule_info.h5
    │   │   │   ├── raw_feature_bc_matrix.h5
    │   │   │   └── web_summary.html
    │   │   ├── WT7_MMT_cellranger
    │   │   │   ├── filtered_feature_bc_matrix.h5
    │   │   │   ├── metrics_summary.csv
    │   │   │   ├── molecule_info.h5
    │   │   │   ├── raw_feature_bc_matrix.h5
    │   │   │   └── web_summary.html
    │   │   ├── WT8_MMT_cellranger
    │   │   │   ├── filtered_feature_bc_matrix.h5
    │   │   │   ├── metrics_summary.csv
    │   │   │   ├── molecule_info.h5
    │   │   │   ├── raw_feature_bc_matrix.h5
    │   │   │   └── web_summary.html
    │   │   └── WT_MMT_cellranger
    │   │       ├── filtered_feature_bc_matrix.h5
    │   │       ├── metrics_summary.csv
    │   │       ├── molecule_info.h5
    │   │       ├── raw_feature_bc_matrix.h5
    │   │       └── web_summary.html
    │   └── master                                            
    │       ├── DGE                                         # Raw differential expression results using the poisson model with speedglm
    │       │   ├── poisson_glm_labelPred_cpc.pdf
    │       │   ├── poisson_glm_labelPred_it_all_0_0_meta.rds
    │       │   ├── poisson_glm_labelPred_it_all_0_0_summary.rds
    │       │   ├── poisson_glm_labelPred_it_all_0_1_meta.rds
    │       │   ├── poisson_glm_labelPred_it_all_0_1_summary.rds
    │       │   ├── poisson_glm_labelPred_it_all_4_0_meta.rds
    │       │   ├── poisson_glm_labelPred_it_all_4_0_summary.rds
    │       │   ├── poisson_glm_labelPred_it_all_6_0_meta.rds
    │       │   ├── poisson_glm_labelPred_it_all_6_0_summary.rds
    │       │   ├── poisson_glm_labelPred_it_all_6_1_meta.rds
    │       │   └── poisson_glm_labelPred_it_all_6_1_summary.rds
    │       ├── clustering_typing                           # results from running clustering, UMAP and findmarkers analyses
    │       │   ├── low_res_clustering_plots_Tcells_harmony.pdf
    │       │   ├── low_res_clustering_plots_harmony.pdf
    │       │   ├── master_seurat_Tcells_low_res_clustering_table_harmony.rds
    │       │   ├── master_seurat_low_res_clustering_table_harmony.rds
    │       │   ├── presto_paired_markers_Sub_res0.05_filtered.rds
    │       │   └── presto_paired_markers_Tcells_curated_cluster_filtered.rds
    │       └── seurat_objects                              # Seurat objects after clustering and UMAP computations
    │           ├── master_seurat_harmony_Tcells_low_res_clustered.rds
    │           └── master_seurat_harmony_low_res_clustered.rds
    ├── src
    │   └── cellranger_count.sh                             # cellranger count script. Run on Yale Ruddle cluster with SLURM scheduler
    └── tables                                              # Final DGE results tables
        ├── Treg_poisson_glm_DGE_results.csv
        ├── central_Tconv_poisson_glm_DGE_results.csv
        ├── effector_Tconv_poisson_glm_DGE_results.csv
        └── resting_Tconv_poisson_glm_DGE_results.csv
```


