#' ---
#' title: "EAE xp: clustering and typing"
#' author: "PPA"
#' date: "`r Sys.Date()`"
#' output:
#'   pdf_document: 
#'     number_sections: yes
#'     toc: yes
#'     toc_depth: 4
#'   html_document:
#'     df_print: paged
#' geometry: margin=0.5in
#' header-includes:
#'    - \usepackage{subfig}
#'    - \usepackage{float}
#' ---
#' 
#' 
## ----setup, include=FALSE--------------------------------------------------------------

##########################
## Packages and options ##
##########################

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
library(forcats)
library(kableExtra)
library(sp)
# library(FactoMineR)
# library(factoextra)
library(wesanderson)

knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.pos = 'H')
knitr::opts_chunk$set(dev = 'png', dpi = 300)

#set root directory to be one level up of the current directory
#knitr::opts_knit$set(root.dir = '.')

#Caching option when tweaking the knitting
knitr::opts_chunk$set(cache = T)

WRITE_FILES <- FALSE


##################
## ggplot theme ##
##################

theme_mrl <- function(x = 1, ...) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(size = 12*x,face = "bold", angle = 45, vjust = 0.9, hjust=0.9, color = "black"),
      axis.text.y = element_text(size = 12*x,face = "bold", color = "black"),
      axis.title.x = element_text(size = 12*x,face = "bold"),
      axis.title.y = element_text(size = 12*x,face = "bold"),
      strip.background = element_rect(fill="gray20", colour="gray20", linetype="solid"),
                                      strip.text = element_text(size=12*x, colour="white", face="bold"),
                                      legend.title = element_text(size=14*x, face = "bold"),
                                      legend.text = element_text(size=12*x, color="black", face="bold"),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
      plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
      plot.caption = element_text(hjust = 0, face = "italic")
                                      ) +
    theme(...)
}


opts <- options()  # save old options

options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

palette <- c(pal_rickandmorty()(12),pal_jco()(10), pal_futurama()(11), pal_tron()(7))




#' 
#' 
#' 
#' 
## ---- fig.width=10---------------------------------------------------------------------
# function to overlay values onto the UMAP embeddings
hex_plots <- function(data,dim1,dim2, color, facets = "", CI = 0.05, fun = "median") {

  message("calc CI")
  # calculate the quantilles for CI limits onthe color scale.
  if(CI) {
    lims <- quantile(pull(data,color), probs = c(CI/2,1 - CI/2))
  } else {
    lims <- c(NA,NA)
  }
  
  message("compute plot")
  ggplot(data, aes_string(dim1, dim2)) +
  stat_summary_hex(bins = 100,
                   fun = fun,
                   aes(z= !!ensym(color), color = ..value..)) +
    facet_wrap(as.formula(paste(". ~", facets)), scales = "free") +
    ggtitle(color) +
    theme_mrl() +
    guides(fill = guide_none()) +
    scale_fill_viridis(limits = lims, oob = squish) +
    scale_color_viridis(limits = lims, oob = squish)
}



ridges_plots <- function(data,dim1,dim2) {
  
  ggplot(data, 
       aes_string(dim1, 
                  dim2)) +
  stat_density_ridges(aes(fill = factor(stat(quantile))),
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = c(0.01,0.05, 0.1,0.9,0.95,0.99), quantile_lines = TRUE,
    color = NA) +
  facet_grid(curated_types~., scales = "free_y", space = "free_y") +
  scale_fill_manual(labels = c("0.01","0.05", "0.1","0.1-0.9","0.9","0.95","0.99"),
                    values = c("black","grey80","black","grey80","black","grey80","black"))
}








#' 
#' 
#' 
#' # clustering and embeddings
#' 
#' # Basic QC metrics
#' 
#' 
## ---- fig.width=12, fig.asp=1.3--------------------------------------------------------

meta_clustering <- readRDS("results/master/clustering_typing/master_seurat_low_res_clustering_table_harmony.rds")


p0 <- ggplot(meta_clustering, aes(Sub_res0.05)) +
  geom_bar() +
  geom_label(stat='count', aes(label=..count..), vjust=1.1) +
  scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p1 <- ggplot(meta_clustering, aes(Sub_res0.05, fill = paste(Geno,sample), color = Geno)) +
  geom_bar(position = "fill", size = 2) +
  scale_color_manual(values = c("red","black")) +
  #scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p2 <- ggplot(meta_clustering, aes(Sub_res0.05, nCount_RNA)) +
  geom_violin() +
  scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p3 <- ggplot(meta_clustering, aes(Sub_res0.05, nFeature_RNA)) +
  geom_violin() +
  scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p4 <- ggplot(meta_clustering, aes(Sub_res0.05, percent.mt)) +
  geom_violin() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p5 <- ggplot(meta_clustering, aes(Sub_res0.05, percent.ribo)) +
  geom_violin() +
  theme_mrl() +
  theme(
        legend.position = "none")


plot_grid(p0,p1,p2,p3,p4,p5, align = "vh", axis = "blrt", ncol = 1)




#' 
#' 
#' 
#' ## cluster split by sample
#' 
#' 
## ----fig.width=12, fig.asp=1-----------------------------------------------------------
p0 <- ggplot(meta_clustering, aes(paste(Xp,sample))) +
  geom_bar() +
  #scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p1 <- ggplot(meta_clustering, aes(paste(Xp,sample), fill = Sub_res0.05, color = Geno)) +
  geom_bar(position = "fill", size = 2) +
  scale_color_manual(values = c("red","black")) +
  #scale_y_log10() +
  theme_mrl()


plot_grid(p0,p1, align = "vh", axis = "blrt", ncol = 1)

#' 
#' 
#' ## Antigen receptor chain expression
#' 
## ---- fig.width=12---------------------------------------------------------------------



p1 <- ggplot(meta_clustering, aes(Sub_res0.05, percent.TCRab+ 0.01)) +
  geom_violin(scale = "width",
              draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_log10()


p2 <- ggplot(meta_clustering, aes(Sub_res0.05, percent.TCRgd+ 0.01)) +
  geom_violin(scale = "width",
              draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_log10()

p3 <- ggplot(meta_clustering, aes(Sub_res0.05, percent.BCR+ 0.01)) +
  geom_violin(scale = "width",
              draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_log10()


p4 <- plot_grid(p1,p2,p3, align = "vh", axis = "blrt", ncol = 1)


p5 <- ggplot(meta_clustering, aes(percent.TCRab+0.01, percent.TCRgd+ 0.01)) +
  geom_hex(aes(color = ..ncount.., fill = ..ncount..)) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_continuous(trans = "log10") +
  scale_fill_continuous(trans = "log10") +
  facet_wrap(.~Sub_res0.05) +
  geom_abline(color = "red")


plot_grid(p4,p5, ncol = 2)





#' 
#' 
#' # Overlap with UMAP embeddings
#' 
#' 
## ---- fig.width=12---------------------------------------------------------------------



ggplot(meta_clustering, aes_string("UMAPn50xMd05xS5_1","UMAPn50xMd05xS5_2")) +
  geom_hex(aes(fill = SNNk20_algo1_res0.1, alpha = log(..ndensity..)),
           bins = 200) +
  scale_alpha(range = c(0.3,1)) +
  theme_mrl()


ggplot(meta_clustering, aes_string("UMAPn50xMd05xS5_1","UMAPn50xMd05xS5_2")) +
  geom_hex(aes(fill = Sub_res0.05, alpha = log(..ndensity..)),
           bins = 200) +
  scale_alpha(range = c(0.3,1)) +
  theme_mrl()



ggplot(meta_clustering, aes_string("UMAPn50xMd05xS5_1","UMAPn50xMd05xS5_2")) +
  geom_hex(aes(fill = Sub_res0.05, alpha = log(..ndensity..)),
           bins = 70) +
  scale_alpha(range = c(0.3,1)) +
  theme_mrl() +
  facet_wrap(.~Sub_res0.05)




meta_clustering %>%
  gather("feat", "val",  -!matches("percent|_RNA$")) %>%
  mutate(val = ifelse(grepl("_RNA$",feat), log10(val), val)) %>%
  group_by(feat) %>%
  nest() %>%
  mutate(plots = purrr::map2(data,feat, ~ hex_plots(.x,
          "UMAPn50xMd05xS5_1",
          "UMAPn50xMd05xS5_2",
          "val",
          ".",
          CI = 0.05) + 
            ggtitle(.y) +
            theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  panel.background = element_rect(fill = 'grey90', color = NA),
                  panel.border = element_rect(colour = "black", fill=NA),
                  strip.text = element_blank(),
                  legend.position = c(1,0),
                  legend.justification = c(1,0),
                  legend.title = element_blank(),
                  legend.direction = "horizontal",
                  legend.text = element_text(angle = 45, hjust = 1, vjust = 1.2))
          )
         ) %>%
  pull(plots) %>%
  plot_grid(plotlist = .)







#' 
#' 
#' 
#' 
#' # Gene cluster typing
#' 
#' 
#' ## Markers
#' 
## --------------------------------------------------------------------------------------

# presto_paired <- readRDS("results/clustering_typing/presto_paired_markers_Sub_res0.05.rds")
# 
# 
# presto_paired <- do.call(rbind, presto_paired) %>%
#   dplyr::filter(padj < 0.01 & logFC > 0 & auc > 0.6)
  
presto_paired <- readRDS("results/master/clustering_typing/presto_paired_markers_Sub_res0.05_filtered.rds")

presto_paired %>%
  group_by(group, feature) %>%
  summarize(count = n(),
            logFC = mean(logFC),
            mlog10pdj = mean(-log10(padj+ 10^-250))) %>%
  dplyr::filter(count > 3) %>% 
  View()



#' 
#' 
#' 
#' from having a look at marker genes and signatures:
#' 
#' * 0_0: CD4 T cells: CD4,CD40lg,cd28
#' * 0_1: CD8  T cells: Cd8a, Cd8b1, Fasl, Gzmb, Gzmk, Klrc1, Klrd1, Nkg7, ccl5,ccl4
#' * 1_0: mixture of cell innate like and Tgd cells. RORC, Cxcr6,Cd3g,Ikzf3,Il17a,Il1r1, il23r, TCR gamma chains.
#' * 2_0: Treg: ikzf2, Foxp3, ctla4, tigit il2ra, Icos, Ass1 (enzyme up in Treg see VJ paper on metabolic control)
#' * 2_1: very small cluster split off of Tregs, with some lymphoid markers present. It carries stress/mitosis features.
#' * 3_0: bcl2, ccl5, Cd7, Nkg7, Tcf7, Sell, lots of ribo genes. Seems transcriptionally very active but not a doublet or bad quality cell. It's a big cluster too. FAU is a potential Tfh marker but no Cd4 marker. Maybe Cd8 cells given it's close to the main CD8 cluster... more central mem given there's no Gzm. Being less stringent on auc threshold we see Cd8 genes and also cd226, Cd160, Pdcd4 (supposedly progenitor marker).
#' * 4_0 and 4_1: B cells: Cd19, Ms4a1, Cd74, Cd79a, CD79b, CD37, MHCgenes, Mzb1 => marginal zone, Bank1
#' * 5_0: Birc5 survival genes, bunch of histones genes, Hmgb's genes: damaged cells?
#' * 5_1: similar less promitent
#' * 6_0: myeloid, CD14, CD74, CD68, CD300c, Tyrobp, Lyz2, Apoe
#' * 6_1: myeloid, ccl4, Tyrobp, CD164, CD180,CD300c, CD47, CD68 MHCgenes
#' 
#' 
#' ## Frequencies of clusters across samples / genotypes
#' 
## ---- fig.width=7----------------------------------------------------------------------

clusters_freqs <- meta_clustering %>%
  #dplyr::filter(!grepl("3_1|^[4-6]",Sub_res0.05)) %>%
  group_by(sample, Geno, Xp, fine) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) 


clusters_freqs %>%
  ggplot(aes(paste(Geno, sample), freq, fill = Sub_res0.05)) +
  geom_bar(stat = "identity") +
  theme_mrl()


clusters_freqs %>%
  ggplot(aes(Geno, freq)) +
  geom_jitter(aes(shape = Xp, color = grepl("WT2|KI8",sample)),
              width = 0.2) +
  stat_summary() +
  facet_wrap(.~Sub_res0.05, scales = "free")




#' 
#' 
#' 
#' 
#' # Subsetting to T cells and co and rerun PCA/harmony/clustering/embeddings
#' 
#' 
#' # Basic QC metrics
#' 
#' 
## ---- fig.width=15, fig.asp=1.3--------------------------------------------------------

meta_clustering <- readRDS("results/master/clustering_typing/master_seurat_Tcells_low_res_clustering_table_harmony.rds")

meta_clustering <- meta_clustering %>%
  mutate(SNNk20_algo1_res0.3 = as.character(SNNk20_algo1_res0.3)) %>%
  mutate(Sub_res0.15 = as.character(Sub_res0.15)) %>%
  mutate(Sub_res0.1 = as.character(Sub_res0.1)) %>%
  mutate(curated_cluster = ifelse(SNNk20_algo1_res0.3 == "6", SNNk20_algo1_res0.3,
                                  ifelse(SNNk20_algo1_res0.3 == "2", Sub_res0.15, Sub_res0.1)))

p0 <- ggplot(meta_clustering, aes(curated_cluster)) +
  geom_bar() +
  geom_label(stat='count', aes(label=..count..), vjust=1.1) +
  scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p1 <- ggplot(meta_clustering, aes(curated_cluster, fill = paste(Geno,sample), color = Geno)) +
  geom_bar(position = "fill", size = 2) +
  scale_color_manual(values = c("red","black")) +
  #scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p2 <- ggplot(meta_clustering, aes(curated_cluster, nCount_RNA)) +
  geom_violin() +
  scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p3 <- ggplot(meta_clustering, aes(curated_cluster, nFeature_RNA)) +
  geom_violin() +
  scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p4 <- ggplot(meta_clustering, aes(curated_cluster, percent.mt)) +
  geom_violin() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p5 <- ggplot(meta_clustering, aes(curated_cluster, percent.ribo)) +
  geom_violin() +
  theme_mrl() +
  theme(
        legend.position = "none")


plot_grid(p0,p1,p2,p3,p4,p5, align = "vh", axis = "blrt", ncol = 1)




#' 
#' 
#' 
#' 
#' 
#' ## cluster split by sample
#' 
#' 
## ----fig.width=15, fig.asp=1-----------------------------------------------------------
p0 <- ggplot(meta_clustering, aes(paste(Xp,sample))) +
  geom_bar() +
  #scale_y_log10() +
  theme_mrl() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p1 <- ggplot(meta_clustering, aes(paste(Xp,sample), fill = curated_cluster, color = Geno)) +
  geom_bar(position = "fill", size = 2) +
  scale_color_manual(values = c("red","black")) +
  #scale_y_log10() +
  theme_mrl()


plot_grid(p0,p1, align = "vh", axis = "blrt", ncol = 1)

#' 
#' 
#' 
#' ## Antigen receptor chain expression
#' 
## ---- fig.width=15---------------------------------------------------------------------



p1 <- ggplot(meta_clustering, aes(curated_cluster, percent.TCRab+ 0.01)) +
  geom_violin(scale = "width",
              draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_log10()


p2 <- ggplot(meta_clustering, aes(curated_cluster, percent.TCRgd+ 0.01)) +
  geom_violin(scale = "width",
              draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_log10()

p3 <- ggplot(meta_clustering, aes(curated_cluster, percent.BCR+ 0.01)) +
  geom_violin(scale = "width",
              draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_log10()


p4 <- plot_grid(p1,p2,p3, align = "vh", axis = "blrt", ncol = 1)


p5 <- ggplot(meta_clustering, aes(percent.TCRab+0.01, percent.TCRgd+ 0.01)) +
  geom_hex(aes(color = ..ncount.., fill = ..ncount..)) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_continuous(trans = "log10") +
  scale_fill_continuous(trans = "log10") +
  facet_wrap(.~curated_cluster) +
  geom_abline(color = "red")


plot_grid(p4,p5, ncol = 2)





#' 
#' 
#' ## Overlap with UMAP embeddings
#' 
#' 
## ---- fig.width=15---------------------------------------------------------------------

umap_setting <- paste0("UMAPn50xMd05xS10_",
                       c(1,2))

meta_clustering <- meta_clustering %>%
  mutate(SNNk20_algo1_res0.3 = as.character(SNNk20_algo1_res0.3)) %>%
  mutate(Sub_res0.15 = as.character(Sub_res0.15)) %>%
  mutate(Sub_res0.1 = as.character(Sub_res0.1)) %>%
  mutate(curated_cluster = ifelse(SNNk20_algo1_res0.3 %in% c("0","1"), Sub_res0.15, Sub_res0.1))

ggplot(meta_clustering, aes_string(umap_setting[1], umap_setting[2])) +
  geom_hex(bins = 200, aes(color = ..count..))


ggplot(meta_clustering, aes_string(umap_setting[1], umap_setting[2])) +
  geom_hex(aes(fill = curated_cluster, alpha = log(..ndensity..)),
           bins = 200) +
  geom_density2d(aes(group = curated_cluster), bins =4, contour_var = "ndensity", color = "black", size = 2) +
  geom_density2d(aes(color = curated_cluster), bins =4, contour_var = "ndensity") +
  scale_alpha(range = c(0.2,1)) +
  theme_mrl(panel.background = element_rect(fill = 'grey90', color = NA)) +
  #theme() +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette)



ggplot(meta_clustering, aes_string(umap_setting[1], umap_setting[2])) +
  geom_hex(aes(fill = curated_cluster, alpha = log(..ndensity..)),
           bins = 70) +
  geom_density2d(aes(group = curated_cluster), bins =4, contour_var = "ndensity", color = "black", size = 2) +
  geom_density2d(aes(color = curated_cluster), bins =4, contour_var = "ndensity") +
  scale_alpha(range = c(0.3,1)) +
  theme_mrl() +
  facet_wrap(.~SNNk20_algo1_res0.3) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette)




meta_clustering %>%
  gather("feat", "val",  -!matches("percent|_RNA$")) %>%
  mutate(val = ifelse(grepl("_RNA$",feat), log10(val), val)) %>%
  group_by(feat) %>%
  nest() %>%
  mutate(plots = purrr::map2(data,feat, ~ hex_plots(.x,
          umap_setting[1], 
          umap_setting[2],
          "val",
          ".",
          CI = 0.05) + 
            ggtitle(.y) +
            theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  panel.background = element_rect(fill = 'grey90', color = NA),
                  panel.border = element_rect(colour = "black", fill=NA),
                  strip.text = element_blank(),
                  legend.position = c(1,0),
                  legend.justification = c(1,0),
                  legend.title = element_blank(),
                  legend.direction = "horizontal",
                  legend.text = element_text(angle = 45, hjust = 1, vjust = 1.2))
          )
         ) %>%
  pull(plots) %>%
  plot_grid(plotlist = .)






#' 
#' 
#' 
#' # cluster and cell types representations for paper:
#' 
#' 
#' 
#' 
## --------------------------------------------------------------------------------------
umap_setting <- paste0("UMAPn50xMd05xS10_",
                       c(1,2))




ggplot(meta_clustering, aes_string(umap_setting[1], umap_setting[2])) +
  geom_hex(aes(fill = curated_cluster, alpha = log(..ndensity..)),
           bins = 200) +
  geom_density2d(aes(group = curated_cluster), bins =4, contour_var = "ndensity", color = "black", size = 2) +
  geom_density2d(aes(color = curated_cluster), bins =4, contour_var = "ndensity") +
  scale_alpha(range = c(0.2,1)) +
  theme_mrl(panel.background = element_rect(fill = 'grey90', color = NA)) +
  #theme() +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette)




#' 
#' 
#' 
#' 
#' 
#' 
#' ## Marker genes for each cluster
#' 
#' 
## --------------------------------------------------------------------------------------

# presto_paired <- readRDS("results/master/clustering_typing/presto_paired_markers_Tcells_curated_cluster.rds")
# 
# 
# presto_paired <- do.call(rbind, presto_paired) %>%
#   dplyr::filter(padj < 0.01 & logFC > 0 & auc > 0.6)

presto_paired <- readRDS("results/master/clustering_typing/presto_paired_markers_Tcells_curated_cluster_filtered.rds")
  


presto_paired %>%
  group_by(group, feature) %>%
  summarize(count = n(),
            logFC = mean(logFC),
            mlog10pdj = mean(-log10(padj+ 10^-250))) %>%
  dplyr::filter(count > 3) %>% 
  View()

marker_genes_to_plot <- "Trdv2-2, Nr4a1, Ctla2a, Ly6c2, Lgals1, Xcl1, Ncr1, Ccr2, Cd2, Cd4, C3g, Cd28, Cd69, Hif1a, Icos, Ifng, Il2r, Ly6a, nfkbia, nkg7, Cd40lg, Tigit, Gzm, Cd27, klrc1, nkg7, cd3e, Cd3e, Ccl4, Cd8a, Cd8b1, Cxcr3, Eomes, Tbx21, Trdv4, CD164, CD247, IL17a, Pdcd1, Lmo4, Ppp1r3b, Icos, Cd40l, Cxcr6, Zbtb16,Itgb7, Ly6a, Il7r, Cd44, Icos, Tmem176b, Cd4, Cd5, Foxp3, Il10ra, Ass1, Ctla4, Il2ra, Icos, Ikzf2, Sell, Cd40lg, Cd4$, Icos, Tcf7, Cd28, Il7r, S1pr1, Cd8a, Fasl, Bcl2, Cd7, Klrd1, Ccl5, Nkg7, Fcer1g, Sell, Tcf7, Cd8b, Ccr7, Klf2, Tcf7, Sell, CD8a, Lef1, Tigit, Cd4, Il1r2, Cd40lg, Cd6$, Cd5$, Ccl1, Lag3, Tgfb1, Ccr8, Junb, Jund, csf2, ifng, rorc, Car2, Cd7$, Cd9$, Cll5, Gzmb, klrb1c, Il12rb2, Tbx21" %>% gsub(", |,","|^", .) %>% gsub(",| ", "_____", .)


marker_genes_to_plot <- presto_paired %>%
  dplyr::filter(grepl(marker_genes_to_plot, feature, ignore.case = T)) %>%
  pull(feature) %>%
  unique() %>%
  str_sort()

#' 
#' ## manual annotations
#' 
#' * 0_0: CD4 T cells: Cd4, Cd5, Cd28, Cd69, Cd40lg
#' * 0_1: CD4 T cells: Cd40lg, Cd4, Icos, Tcf7, Cd28, Il7r, S1pr1
#' * 1_0: CD8 T cells: CD8a, CD8b1, all Gzm, Prf1, klrc1, nkg7, cd3e, Cd3e, Ccl4, Cd8a, Cd8b1, Cxcr3, Hopx
#' * 1_1: CD8 T cells: Cd8a, Cd8b1, Ccl5, Itga4, Il10ra, Tigit, Tox, Bcl2, Cd27
#' * 1_2: CD8 T cells: lots of ISGs, Eomes, Il10ra, Ccl5, Gzmb
#' * 2_0: NKT cells: expresses mostly TCRab but a bit of TCrgd too. Lmo4, Ppp1r3b, Icos, Cd40l, Cxcr6, Zbtb16 encodes PLZF master TF of NKT.
#' * 2_1: Trdv2-2 gd T cells: TCR gd segments, Itgb7, Ly6a, Cd82, Il7r, Cd44, Icos, Rorc, Il17a, Tmem176b (marker of Th17 phenotype)
#' 
#' * 3_0: unclear cytotoxic T cells: Bcl2, Cd7, Klrd1 and lots of Klr's, Ccl5, Fcer1g, Sell, Tcf7, Xcl1, Il2rb,  lots of ribo genes expression (maybe prolif).
#' * 3_1: unclear cytotoxic T cells, MAIT or NK: Car2, Cd7, Cd9, Klr's, Ccl5, Gzmb, klrb1c (CD161), Il12rb, Xcl1
#' * 3_2: CD8 T cells: Cd8a, Cd8b, Ccr7, Klf2, Tcf7, Sell, Lef1.
#' 
#' * 4_0: Treg: Cd4, Cd5, Foxp3, Il10ra, Ass1, Ctla4, Il2ra, Icos, Ikzf2, Sell (resting?)
#' 
#' * 5_0: Trdv4 gd T cells: TCR gd segments, CD164, CD247,IL17a, Pdcd1, Rorc, Cd3g
#' 
#' 
#' * 6_0: CD4 T cells: Cd4, Il1r2, Cd40lg, Cd6, Cd5, Ccl1, Lag3, Tgfb1, Ccr8, Junb, Jund, csf2 (encodes gm-csf), ifng, Ccr8, Hif1a
#' 
#' * 6_1: stress/DNA damage/IFN T cells
#' 
#' 
#' 
#' 
#' # cluster and cell types representations for paper:
#' 
#' 
#' 
#' 
## --------------------------------------------------------------------------------------
umap_setting <- paste0("UMAPn50xMd05xS10_",
                       c(1,2))


labels <- data.frame(
  `0_0` = c("central_memory_Tconv", "Tconv"),
  `0_1` = c("resting_Tconv", "Tconv"),
  `1_0` = c("central_memory_CTL", "CTL"),
  `1_1` = c("effector_memory_CTL", "CTL"),
  `1_2` = c("IFNs_stimulated_CTL", "CTL"),
  `2_0` = c("iNKT", "iNKT"),
  `2_1` = c("Trdv2-2_gdT", "gdT"),
  `3_0` = c("NKT-like", "NKT-like"),
  `3_1` = c("NK", "NK"),
  `3_2` = c("resting_CTL", "CTL"),
  `4_0` = c("Treg", "Treg"),
  `5_0` = c("Trdv4_gdT", "gdT"),
  `6_0` = c("effector_memory_Tconv", "Tconv"),
  `6_1` = c("stressed_T_cells", "stressed_T_cells"),
  row.names = c("fine", "main")
) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("curated_cluster") %>%
  mutate(curated_cluster = gsub("^X","", curated_cluster))



meta_clustering %>%
  left_join(labels) %>%
ggplot(aes_string(umap_setting[1], umap_setting[2])) +
  geom_point(aes(fill = main), shape = 21, alpha = 0.2, size = 1.5, stroke = 0) +
  theme(panel.background = element_rect(fill = 'grey90', color = NA),
        panel.border = element_rect(colour = "black", fill=NA, linetype = "solid", size = 1),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        aspect.ratio = 0.7) +
  scale_fill_igv() +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 5, alpha = 1)))


ggsave("figures/master/clusters/main.png", height = 5, width = 7)
ggsave("figures/master/clusters/main_legend.pdf", get_legend(last_plot()),height = 4, width = 3)


meta_clustering %>%
  left_join(labels) %>%
  mutate(fine = gsub("_"," ", fine)) %>%
ggplot(aes_string(umap_setting[1], umap_setting[2])) +
  geom_point(aes(fill = fct_reorder(fine, rank(main))), shape = 21, alpha = 0.2, size = 1.5, stroke = 0) +
  theme(panel.background = element_rect(fill = 'grey90', color = NA),
        panel.border = element_rect(colour = "black", fill=NA, linetype = "solid", size = 1),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        aspect.ratio = 0.7) +
  scale_fill_manual(values = palette[5:30]) +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 5, alpha = 1)))


ggsave("figures/master/clusters/fine.png", height = 5, width = 7)
ggsave("figures/master/clusters/fine_legend.pdf", get_legend(last_plot()),height = 4, width = 3)



#' 
#' 
#' 
#' 
#' ## Frequencies of clusters across samples / genotypes
#' 
## ---- fig.width=15---------------------------------------------------------------------

clusters_freqs <- meta_clustering %>%
  #dplyr::filter(!grepl("3_1|^[4-6]",Sub_res0.05)) %>%
  group_by(sample, Geno, Xp, curated_cluster) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) 


clusters_freqs %>%
  ggplot(aes(paste(Geno, sample), freq, fill = curated_cluster)) +
  geom_bar(stat = "identity") +
  theme_mrl()


new_scale <- scales::trans_new("shift",
                             transform = function(x) {log10(x+1)},
                             inverse = function(x) {10^x - 1})


clusters_freqs %>%
  left_join(labels) %>%
  mutate(fine = gsub("_"," ", fine)) %>%
  ggplot(aes(fct_rev(Geno), freq*100, fill = fct_rev(Geno))) +
  stat_summary(aes(width = 0.5),fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_summary(fun.data = mean_se, geom = "bar", color = "black") +
  scale_fill_manual(values = c("grey","#476749")) +
  ggsignif::geom_signif(comparisons = list(c("KI", "WT")), map_signif_level = T) +
  facet_wrap(.~fct_reorder(fine, rank(main)), nrow = 1, labeller = label_wrap_gen(10), strip.position = "bottom") +
  scale_y_continuous(trans = new_scale, 
                     expand = expand_scale(c(0,0.1)),
                     breaks = c(0.3,1,3,10,30)) +
  theme_mrl(0.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(colour = "black", angle = 45),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0, "lines")) +
  labs(y = "%", x = NULL) +
  guides(fill = guide_legend(title = "Genotype"))


ggsave("figures/master/clusters/freq.pdf", height = 3, width = 10)



#' 
#' ## generate plots for marker gene expression projected on UMAP
#' 
## --------------------------------------------------------------------------------------


norm_counts <- readRDS("results/master/seurat_objects/master_seurat_harmony_Tcells_low_res_clustered.rds")[["RNA"]]@data

norm_counts <- norm_counts[rownames(norm_counts) %in% marker_genes_to_plot,] %>% #marker_genes_to_plot
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename_with(~paste0("gene_",.))


meta_clustering %>%
  dplyr::select(cell.id,UMAPn50xMd05xS10_1,UMAPn50xMd05xS10_2, curated_cluster) %>%
  left_join(norm_counts %>% tibble::rownames_to_column("cell.id")) %>%
  gather("feat", "val",  -!matches("^gene_")) %>%
  group_by(feat) %>%
  nest() %>%
  mutate(plots = purrr::map2(data,feat, ~ ggsave(paste0("figures/master/marker_genes/",.y,"_mean.png"),
                                                 hex_plots(.x,
                                                           "UMAPn50xMd05xS10_1",
                                                           "UMAPn50xMd05xS10_2",
                                                           "val",
                                                           ".",
                                                           CI = 0.05,
                                                           fun = "mean") + 
                                                   ggtitle(.y) +
            theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  panel.grid = element_blank(),
                  panel.background = element_rect(fill = 'grey90', color = NA),
                  panel.border = element_rect(colour = "black", fill=NA),
                  strip.text = element_blank(),
                  aspect.ratio = 0.7,
                  legend.position = c(1,0),
                  legend.justification = c(1,0),
                  legend.title = element_blank(),
                  legend.direction = "horizontal",
                  legend.text = element_text(angle = 45, hjust = 1, vjust = 1.2)))
          )
         )
  


#' 
