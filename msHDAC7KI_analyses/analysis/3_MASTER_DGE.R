#' ---
#' title: "EAE xp: DGE"
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
library(Matrix)
library(DESeq2)
library(edgeR)
library(ashr)
library(vsn)

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



# import gene annotations

gene_annotations <- readRDS("data/gene_annotations_positions.rds")



#' 
#' 
#' 
#' 
#' # poisson glm
#' 
#' 
## ---- fig.width= 10--------------------------------------------------------------------


files <- list.files("results/master/DGE/", "*summary.rds", full.names = F)

poissonglm_results <- lapply(files, function(file) {
  
  
  
  res <- readRDS(file.path("results/master/DGE/",file)) %>%
    mutate(iteration = gsub("poisson_glm_labelPred_it_|_summary.rds","",file))
  

  
})

poissonglm_results <- do.call(rbind, poissonglm_results) %>%
  separate(iteration, c("iteration", "curated_cluster"), sep = "_", extra = "merge") %>%
  dplyr::filter(term == "GenoWT")


# apply shrinkage and fdr
cleaned_poissonglm_results <- poissonglm_results %>%
  group_by(iteration, curated_cluster) %>%
  nest() %>%
  mutate(fdr = purrr::map(data, ~ p.adjust(.x$p.value, method = "fdr"))) %>%
  mutate(ashr_results = purrr::map(data, ~ashr::ash(.x$estimate, .x$std.error, mixcompdist="normal", method="shrink")$result)) %>%
  unnest() %>% 
  rename(sLogFC = PosteriorMean) 





#' 
## --------------------------------------------------------------------------------------
######
# add gene expression levels
files <- list.files("results/master/DGE/", "*meta.rds", full.names = F)

poissonglm_meta <- lapply(files, function(file) {
  
  
  
  res <- readRDS(file.path("results/master/DGE/",file)) %>%
    mutate(iteration = gsub("poisson_glm_labelPred_it_|_meta.rds","",file))
  

  
})



poissonglm_meta <- do.call(rbind, poissonglm_meta) %>%
  separate(iteration, c("iteration", "curated_cluster"), sep = "_", extra = "merge")

cleaned_poissonglm_results <- cleaned_poissonglm_results %>%
  left_join(poissonglm_meta) %>%
  group_by(iteration, curated_cluster) %>%
  nest()


#' 
#' 
#' ## Generate volcano, ma, p hist plots
#' 
## ---- fig.width= 20--------------------------------------------------------------------


# save volcano plots
cleaned_poissonglm_results %>%
  mutate(plots = purrr::map2(
    data,
    paste(iteration, curated_cluster),
    ~ ggsave(
      paste0("figures/master/DGE/", .y, "volcano_plot.png"),
      ggplot(.x, aes(sLogFC,-log10(p.value + 10 ^ -250))) +
        geom_hex(bins = 100) +
        scale_fill_viridis(trans = "log10") +
        geom_text_repel(
          data = dplyr::filter(.,!grepl("^Tra|^Trb|^Trd|^Trg|^Gm|Rik$|^Rps|^Rpl",gene_name)) %>% dplyr::filter(rank(-abs(sLogFC)) < 50 | rank(fdr) < 50),
          aes(label = gene_name),
          color = "red",
          size = 6
        ) +
        labs(title = .y,
             caption = sum(.x$fdr < 0.01))
    )
  ))


# cleaned_poissonglm_results <- poissonglm_results %>%
#   group_by(iteration, curated_cluster) %>%
#   nest() %>%
#   mutate(fdr = purrr::map(data, ~ p.adjust(.x$p.value, method = "fdr"))) %>%
#   mutate(ashr_results = purrr::map(data, ~ashr::ash(.x$estimate, .x$std.error, mixcompdist="normal", method="fdr")$result)) %>%
#   unnest() %>% 
#   rename(sLogFC = PosteriorMean) %>%
#   group_by(iteration, curated_cluster) %>%
#   nest()

# volcano plots with ashr PP genes
cleaned_poissonglm_results %>%
  mutate(plots = purrr::map2(
    data,
    paste(iteration, curated_cluster),
    ~ ggsave(
      paste0("figures/master/DGE/", .y, "volcano_plot_ashr_fdr.png"),
      ggplot(.x, aes(sLogFC,-log10(p.value + 10 ^ -250))) +
        geom_hex(bins = 100) +
        scale_fill_viridis(trans = "log10") +
        geom_point(
          data = dplyr::filter(
            .,
            !grepl("^Tra|^Trb|^Trd|^Trg|^Gm|Rik$|^Rps|^Rpl", gene_name)
          ) %>% dplyr::filter(PositiveProb > 0.95 | NegativeProb > 0.95),
          color = "red"
        ) +
        geom_text_repel(
          data = dplyr::filter(
            .,
            !grepl("^Tra|^Trb|^Trd|^Trg|^Gm|Rik$|^Rps|^Rpl", gene_name)
          ) %>% dplyr::filter(rank(-abs(sLogFC)) < 50 | rank(fdr) < 50),
          aes(label = gene_name),
          color = "red",
          size = 6,
          force = 5
        ) +
        labs(title = .y,
             caption = sum(
               .x$fdr < 0.1 & (.x$NegativeProb > 0.95 | .x$PositiveProb > 0.95)
             )) +
        facet_grid(p.value > quantile(.$p.value, 0.01) ~ ., scales = "free_y")
    )
  ))



# MA plots with ashr PP genes
cleaned_poissonglm_results %>%
  mutate(plots = purrr::map2(
    data,
    paste(iteration, curated_cluster),
    ~ ggsave(
      paste0("figures/master/DGE/", .y, "ma_plot_ashr_fdr.png"),
      ggplot(.x, aes(log10(cpc),sLogFC)) +
        geom_hex(bins = 70) +
        scale_fill_viridis(trans = "log10") +
        geom_point(
          data = dplyr::filter(
            .,
            !grepl("^Tra|^Trb|^Trd|^Trg|^Gm|Rik$|^Rps|^Rpl", gene_name)
          ) %>% dplyr::filter(PositiveProb > 0.95 | NegativeProb > 0.95),
          color = "red"
        ) +
        geom_label_repel(
          data = dplyr::filter(
            .,
            !grepl("^Tra|^Trb|^Trd|^Trg|^Gm|Rik$|^Rps|^Rpl", gene_name)
          ) %>% dplyr::filter(rank(-abs(sLogFC)) < 50 | rank(fdr) < 50),
          aes(label = gene_name),
          color = "red",
          size = 6,
          force = 5
        ) +
        labs(title = .y,
             caption = sum(
               .x$fdr < 0.1 & (.x$NegativeProb > 0.95 | .x$PositiveProb > 0.95)
             )) +
        scale_y_continuous(limits = c(-1,1), oob = squish) +
        theme(plot.background = element_rect(fill="white")),
        
    )
  ))

# histogram p values for each iteration 
 cleaned_poissonglm_results %>%
  mutate(plots = purrr::map2(
    data,
    paste(iteration, curated_cluster),
    ~ ggsave(
      paste0("figures/master/DGE/", .y, "p_val_histogram.png"),
      ggplot(.x, aes(p.value)) +
        geom_histogram() +
        labs(title = .y,
             caption = sum(.x$fdr < 0.1))
    )
  ))

 
 
 

#' 
#' # shrinkage effects
#' 
#' 
## ---- fig.width= 10--------------------------------------------------------------------


# compare parameters between iterations

beta_assemble <- function(dat, new_name) {
    dplyr::select(dat, gene_name, curated_cluster,sLogFC) %>% rename(!!new_name := sLogFC)
  }

cleaned_poissonglm_results %>%
  unnest() %>%
  group_by(iteration) %>%
  nest() %>%
  mutate(data = purrr::map2(data, iteration, ~ beta_assemble(.x, .y))) %>%
  pull(data) %>%
  purrr::reduce(full_join, by = c("gene_name", "curated_cluster")) %>%
  dplyr::select(-gene_name) %>%
  group_by(curated_cluster) %>%
  nest() %>%
  mutate(plots = purrr::map2(data, curated_cluster, ~ print(GGally::ggpairs(.x) + 
                                                              ggtitle(.y))))
  


several_iterations = F

if (several_iterations) {
  library(UpSetR)
  
  add_names <- function(x, name) {
    names(x) <- name
    return(x)
  }
  
  
  cleaned_poissonglm_results %>%
    mutate(DEGs = purrr::map(data, ~ {
      .x %>%
        dplyr::filter(fdr < 0.05) %>%
        pull(gene_name)
      
    })) %>%
    dplyr::select(-data) %>%
    group_by(curated_cluster) %>%
    summarize(DEGs = list(DEGs), iteration =  list(paste(curated_cluster, iteration))) %>%
    mutate(DEGs = purrr::map2(DEGs, iteration, ~ add_names(.x, .y))) %>%
    mutate(plot = purrr::map(DEGs, ~ upset(fromList(.x)))) %>%
    pull(plot) %>%
    lapply(., print)
}

# visualize shrinkage changes for the different clusters

cleaned_poissonglm_results %>%
  unnest() %>%
  dplyr::filter(
    #grepl("X14mKI8", iteration) &
      #curated_cluster %in% c("0_0", "4_0", "6", "3_0") &
      fdr < 0.05
  ) %>%
  {
    ggplot(., aes(estimate, sLogFC)) +
      theme_mrl() +
      geom_hex(bins = 80) +
      facet_wrap(. ~ curated_cluster, scales = "free") +
      scale_fill_viridis(trans = "log10") +
      geom_abline(color = "red") +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_text(
        data = . %>% group_by(curated_cluster) %>% summarize(count = sum(fdr < 0.05)),
        inherit.aes = F,
        aes(
          x = -Inf,
          y = Inf,
          label = paste(count, "DEGs at fdr < 0.05")
        ),
        hjust = -0.1,
        vjust = 1.5
      )
  }


#' 
#' 
## --------------------------------------------------------------------------------------

cleaned_poissonglm_results %>%
  unnest() %>%
  dplyr::filter(fdr < 0.05 & curated_cluster == "4_0") %>%
  View()






#' 
#' 
#' # Curated plots for manuscript
#' 
#' ## curated plot for Treg cells
#' 
## ---- fig.width=10---------------------------------------------------------------------
genes_to_annotate <- "Ccl1|Ccl4|Ccl5|Il22|Il4|Cd74|Areg|Ifng|Klrg1|Hdac7|Il10ra|Cd74|Cd47|Ccr9|Il17a|Tox2|Maf|Ahr|Ifngr1|Cd8a|Cd8b1|Cd163l1"

genes_to_annotate <- c("Ccl1","Ccl4","Ccl5","Il22","Il4","Areg","Ifng","Klrg1","Hdac7","Il10ra",
                       "Cd74","Cd47","Ccr9","Il17a","Tox2","Maf","Ahr","Ifngr1",
                       "Gzmb", "Hopx", "Tox2", "Tigit", "Havcr2", "Pdcd1", "Cd226","Tnfrsf9", "Lag3", "Ctla4")

cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "4_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  #dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05)
  {ggplot(.,aes(avg_expr, -sLogFC)) +
      geom_point(aes(color = fdr < 0.05)) +
      theme_mrl() +
      scale_y_continuous(limits = c(-1,1), oob = squish) +
      scale_x_log10() +
      scale_color_manual(values = c("grey", "red")) +
      geom_label_repel(data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05),
                       aes(label = gene_name, fill = "red"), alpha = 1, min.segment.length = 0) +
      geom_point(data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05),
                       shape = 21, color = "black", size =2, fill = "red") +
      labs(x = "log10 average gene expression", y = "shrunk Log Fold Change")
    }

cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "4_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  dplyr::filter(fdr < 0.05) %>%
  write_csv("tables/Treg_poisson_glm_DEGs.csv")

cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "4_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::select(
    gene_name,
    p_val = p.value,
    fdr,
    LogFC = estimate,
    sLogFC,
    avg_expr = cpc
  ) %>%
  dplyr::filter(
    !grepl(
      "^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$",
      gene_name,
      ignore.case = T
    )
  ) %>%
  left_join(dplyr::select(gene_annotations, gene_id, gene_name, external_gene_name)) %>%
  #dplyr::filter(fdr < 0.05) %>%
  write_csv("tables/Treg_poisson_glm_DGE_results.csv")







# remove ribosomial, mito and TCR genes for MA plot
cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "4_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::select(
    gene_name,
    p_val = p.value,
    fdr,
    LogFC = estimate,
    sLogFC,
    avg_expr = cpc
  ) %>%
  dplyr::filter(
    !grepl(
      "^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$",
      gene_name,
      ignore.case = T
    )
  ) %>%
  mutate(cat = ifelse(fdr < 0.05 & -sLogFC > 0, "up", 
                      ifelse(fdr < 0.05 & -sLogFC < 0, "dn", "ns"))) %>%
  #dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05)
  {
    ggplot(., aes(avg_expr,-sLogFC)) +
      geom_point(aes(color = cat)) +
      theme_mrl()  +
      stat_density2d(aes(fill = stat(level)), breaks = seq(1,6,0.1), geom = "polygon") +
      scale_fill_gradient(low = "dark grey", high = "black") +
      scale_y_continuous(limits = c(-1, 1), oob = squish) +
      scale_x_log10() +
      scale_color_manual(values = c(ns = "dark grey", up = "red4", dn = "navy")) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        aes(label = gene_name),
        color = "navy",
        alpha = 1,
        force = 30,
        segment.size = 0.1,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(NA, -0.3)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "navy"
      ) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        aes(label = gene_name),
        color = "red4",
        alpha = 1,
        force = 30,
        segment.size = 0.1,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(0.3,NA)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "red4"
      ) +
      labs(x = "log10 average gene expression", y = "shrunk Log Fold Change") +
      theme(axis.text = element_text(color="black"),
            legend.position = "none",
            panel.border = element_rect(color ="black",fill=NA),
            axis.line = element_blank(),
            aspect.ratio = 1)
  }




ggsave("figures/master/DGE/paper_MAplot_Treg.pdf", width = 5, height = 5)


#' 
#' 
#' 
#' 
#' 
#' 
#' 
## ---- fig.width=10---------------------------------------------------------------------
# remove ribosomial, mito and TCR genes for volcano plot
cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "4_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::select(
    gene_name,
    p_val = p.value,
    fdr,
    LogFC = estimate,
    sLogFC,
    avg_expr = cpc
  ) %>%
  dplyr::filter(!grepl("^RP[SL]|^MT-|^TRA[VJ]|^TRB[VDJ]", gene_name, ignore.case = T)) %>%
  #dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05)
  {
    ggplot(., aes(avg_expr,-sLogFC)) +
      geom_point(aes(color = cat)) +
      theme_mrl()  +
      stat_density2d(aes(fill = stat(level)), breaks = seq(1,6,0.1), geom = "polygon") +
      scale_fill_gradient(low = "dark grey", high = "black") +
      scale_y_continuous(limits = c(-1, 1), oob = squish) +
      scale_x_log10() +
      scale_color_manual(values = c(ns = "dark grey", up = "red4", dn = "navy")) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        aes(label = gene_name),
        color = "navy",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(NA, -0.2)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "navy"
      ) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        aes(label = gene_name),
        color = "red4",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(0.2,NA)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "red4"
      ) +
      labs(x = "log10 average gene expression", y = "shrunk Log Fold Change") +
      theme(axis.text = element_text(color="black"),
            legend.position = "none",
            panel.border = element_rect(color ="black",fill=NA),
            axis.line = element_blank(),
            aspect.ratio = 1)
  }




#' 
#' ## curated plot for Th17 cells
#' 
#' ### Cluster 0: main one
#' 
## ---- fig.width=10---------------------------------------------------------------------

genes_to_annotate <- "Ccl1|Ccl4|Ccl5|Il22|Il4|Cd74|Areg|Ifng|Klrg1|Hdac7|Il10ra|Cd74|Cd47|Ccr9|Il17a|Tox2|Maf|Ahr|Ifngr1|Cd8a|Cd8b1|Cd163l1"

genes_to_annotate <- c("Ccl1","Ccl4","Ccl5","Il22","Il4","Areg","Ifng","Klrg1","Hdac7","Sgk1", "S1pr1", "Lgals3", "Il5","Il17f", "Il23r", "Cpt1a",
                       "Il10ra","Cd74","Cd47","Ccr9","Il17a","Tox2","Maf","Ahr","Ifngr1","Il13","Il22","Cxcl3", "Ccr6", "Il2", "Il3",
                       "Rorc", "Tbx21", "Il17a", "Csf2", "Cxcl3", "Il22", "Gzmb", "Casp1", "Stat4", "Il9", "Il10", "Ikzf3",
                       "Lag3", "Cd226", "Pdcd1", "Tigit", "Havcr2", "Pdcd1", "Cd226", "Lrmp", "Ccl3", "Il1rn")

cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "0_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::filter(!grepl("^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$", gene_name, ignore.case = T)) %>%
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  left_join(dplyr::select(gene_annotations, gene_id, gene_name, external_gene_name)) %>%
  #dplyr::filter(fdr < 0.05) %>%
  write_csv("tables/central_Tconv_poisson_glm_DGE_results.csv")



cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "0_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::filter(!grepl("^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$", gene_name, ignore.case = T)) %>%
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  mutate(cat = ifelse(fdr < 0.05 & -sLogFC > 0, "up", 
                      ifelse(fdr < 0.05 & -sLogFC < 0, "dn", "ns"))) %>%
  #dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05)
  {
    ggplot(., aes(avg_expr,-sLogFC)) +
      geom_point(aes(color = cat)) +
      theme_mrl()  +
      stat_density2d(aes(fill = stat(level)), breaks = seq(1.5,6,0.1), geom = "polygon") +
      scale_fill_gradient(low = "dark grey", high = "black") +
      scale_y_continuous(limits = c(-1, 1), oob = squish) +
      scale_x_log10() +
      scale_color_manual(values = c(ns = "dark grey", up = "red4", dn = "navy")) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        aes(label = gene_name),
        color = "navy",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(NA, -0.2)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "navy"
      ) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        aes(label = gene_name),
        color = "red4",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(0.2,NA)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "red4"
      ) +
      labs(x = "log10 average gene expression", y = "shrunk Log Fold Change") +
      theme(axis.text = element_text(color="black"),
            legend.position = "none",
            panel.border = element_rect(color ="black",fill=NA),
            axis.line = element_blank(),
            aspect.ratio = 1)
  }

ggsave("figures/master/DGE/paper_MAplot_central.pdf", width = 5, height = 5)



#' 
#' ### Cluster 6: most effector one
#' 
## ---- fig.width=10---------------------------------------------------------------------

genes_to_annotate <- "Ccl1|Ccl4|Ccl5|Il22|Il4|Cd74|Areg|Ifng|Klrg1|Hdac7|Il10ra|Cd74|Cd47|Ccr9|Il17a|Tox2|Maf|Ahr|Ifngr1|Cd8a|Cd8b1|Cd163l1"

genes_to_annotate <- c("Ccl1","Ccl4","Ccl5","Il22","Il4","Areg","Ifng","Klrg1","Hdac7","Sgk1", "S1pr1", "Lgals3", "Il5","Il17f", "Il23r", "Cpt1a",
                       "Il10ra","Cd74","Cd47","Ccr9","Il17a","Tox2","Maf","Ahr","Ifngr1","Il13","Il22","Cxcl3", "Ccr6", "Il2", "Il3",
                       "Rorc", "Tbx21", "Il17a", "Csf2", "Cxcl3", "Il22", "Gzmb", "Casp1", "Stat4", "Il9", "Il10", "Ikzf3",
                       "Lag3", "Cd226", "Pdcd1", "Tigit", "Havcr2", "Pdcd1", "Cd226", "Lrmp", "Ccl3", "Il1rn")


cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "6_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::filter(!grepl("^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$", gene_name, ignore.case = T)) %>%
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  left_join(dplyr::select(gene_annotations, gene_id, gene_name, external_gene_name)) %>%
  #dplyr::filter(fdr < 0.05) %>%
  write_csv("tables/effector_Tconv_poisson_glm_DGE_results.csv")




cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "6_0") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::filter(!grepl("^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$", gene_name, ignore.case = T)) %>%
  #View()
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  mutate(cat = ifelse(fdr < 0.05 & -sLogFC > 0, "up", 
                      ifelse(fdr < 0.05 & -sLogFC < 0, "dn", "ns"))) %>%
  #dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05)
  {
    ggplot(., aes(avg_expr,-sLogFC)) +
      geom_point(aes(color = cat)) +
      theme_mrl()  +
      stat_density2d(aes(fill = stat(level)), breaks = seq(1,6,0.1), geom = "polygon") +
      scale_fill_gradient(low = "dark grey", high = "black") +
      scale_y_continuous(limits = c(-1, 1), oob = squish) +
      scale_x_log10() +
      scale_color_manual(values = c(ns = "dark grey", up = "red4", dn = "navy")) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        aes(label = gene_name),
        color = "navy",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(NA, -0.1)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "navy"
      ) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        aes(label = gene_name),
        color = "red4",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(0.1,NA)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "red4"
      ) +
      labs(x = "log10 average gene expression", y = "shrunk Log Fold Change") +
      theme(axis.text = element_text(color="black"),
            legend.position = "none",
            panel.border = element_rect(color ="black",fill=NA),
            axis.line = element_blank(),
            aspect.ratio = 1)
  }


ggsave("figures/master/DGE/paper_MAplot_Teff.pdf", width = 5, height = 5)


#' 
#' ### Cluster 0_1: resting memory
#' 
## ---- fig.width=10---------------------------------------------------------------------

genes_to_annotate <- "Ccl1|Ccl4|Ccl5|Il22|Il4|Cd74|Areg|Ifng|Klrg1|Hdac7|Il10ra|Cd74|Cd47|Ccr9|Il17a|Tox2|Maf|Ahr|Ifngr1|Cd8a|Cd8b1|Cd163l1"

genes_to_annotate <- c("Ccl1","Ccl4","Ccl5","Il22","Il4","Areg","Ifng","Klrg1","Hdac7","Sgk1", "S1pr1", "Lgals3", "Il5","Il17f", "Il23r", "Cpt1a",
                       "Il10ra","Cd74","Cd47","Ccr9","Il17a","Tox2","Maf","Ahr","Ifngr1","Il13","Il22","Cxcl3", "Ccr6", "Il2", "Il3",
                       "Rorc", "Tbx21", "Il17a", "Csf2", "Cxcl3", "Il22", "Gzmb", "Casp1", "Stat4", "Il9", "Il10", "Ikzf3",
                       "Lag3", "Cd226", "Pdcd1", "Tigit", "Havcr2", "Pdcd1", "Cd226", "Lrmp", "Ccl3", "Il1rn")


cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "0_1") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::filter(!grepl("^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$", gene_name, ignore.case = T)) %>%
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  left_join(dplyr::select(gene_annotations, gene_id, gene_name, external_gene_name)) %>%
  #dplyr::filter(fdr < 0.05) %>%
  write_csv("tables/resting_Tconv_poisson_glm_DGE_results.csv")



cleaned_poissonglm_results %>%
  dplyr::filter(curated_cluster == "0_1") %>%
  unnest() %>%
  ungroup() %>%
  dplyr::filter(!grepl("^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$", gene_name, ignore.case = T)) %>%
  #View()
  dplyr::select(gene_name, 
                p_val = p.value,
                fdr,
                LogFC = estimate,
                sLogFC,
                avg_expr = cpc) %>%
  mutate(cat = ifelse(fdr < 0.05 & -sLogFC > 0, "up", 
                      ifelse(fdr < 0.05 & -sLogFC < 0, "dn", "ns"))) %>%
  #dplyr::filter(gene_name %in% genes_to_annotate & fdr < 0.05)
  {
    ggplot(., aes(avg_expr,-sLogFC)) +
      geom_point(aes(color = cat)) +
      theme_mrl()  +
      stat_density2d(aes(fill = stat(level)), breaks = seq(1,6,0.1), geom = "polygon") +
      scale_fill_gradient(low = "dark grey", high = "black") +
      scale_y_continuous(limits = c(-1, 1), oob = squish) +
      scale_x_log10() +
      scale_color_manual(values = c(ns = "dark grey", up = "red4", dn = "navy")) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        aes(label = gene_name),
        color = "navy",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(NA, -0.1)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "dn"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "navy"
      ) +
      geom_label_repel(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        aes(label = gene_name),
        color = "red4",
        alpha = 1,
        force = 20,
        min.segment.length = 0,
        segment.color = "black",
        ylim = c(0.1,NA)
      ) +
      geom_point(
        data = . %>% dplyr::filter(gene_name %in% genes_to_annotate & cat == "up"),
        shape = 21,
        color = "black",
        size = 2,
        fill = "red4"
      ) +
      labs(x = "log10 average gene expression", y = "shrunk Log Fold Change") +
      theme(axis.text = element_text(color="black"),
            legend.position = "none",
            panel.border = element_rect(color ="black",fill=NA),
            axis.line = element_blank(),
            aspect.ratio = 1)
  }


ggsave("figures/master/DGE/paper_MAplot_Tresting.pdf", width = 5, height = 5)



