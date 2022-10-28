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
library(ggridges)


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


merged.seurat <- readRDS(file.path(results_dir,"seurat_objects","master_seurat.rds"))





raw_meta <- merged.seurat@meta.data %>%
  tibble::rownames_to_column("cell.id") 
  
  
#################################################
################## QC plots #####################
#################################################

# function to get the highest density points and the mean (output formatted for use in ggplot stat_summary)
high_point <- function(x, high_only = F) {
  den <- density(x)
  high <- data.frame(x = den[["x"]],
                     y = den[["y"]]) %>%
    dplyr::filter(y == max(y)) %>%
    pull(x) %>%
    mean()
  
  if(!high_only) {
    if(mean(x) < high) {
      return(c(y = high,ymin = mean(x), ymax = high))
    } else {
      return(c(y = high,ymin = high, ymax = mean(x)))
    }
  } else {
    return(high)
  }
  
    
  
}
  
  
pdf("doc/master_A02_filter_cells.pdf", 12,9)

#~~~~~~~~~~~~~~~~~~~~~~
#Descriptive plots of computed QC metrics

raw_meta %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(nCount_RNA,nFeature_RNA)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  facet_wrap(.~ paste(Xp,sample)) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks()



raw_meta %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(nCount_RNA,percent.mt)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  facet_wrap(.~ paste(Xp,sample)) +
  scale_x_log10() +
  #scale_y_log10() +
  annotation_logticks(sides = "b") +
  scale_y_continuous(limits = c(NA,10), oob = squish)



raw_meta %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(percent.ribo,percent.mt)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  facet_wrap(.~ paste(Xp,sample)) +
  #scale_x_log10() +
  #scale_y_log10() +
  #annotation_logticks(sides = "b") +
  scale_y_continuous(limits = c(NA,10), oob = squish)


raw_meta %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(nCount_RNA,percent.TCRab)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  facet_wrap(.~ paste(Xp,sample)) +
  scale_x_log10() +
  #scale_y_log10() +
  annotation_logticks(sides = "b")



raw_meta %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(nCount_RNA,percent.TCRgd)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  facet_wrap(.~ paste(Xp,sample)) +
  scale_x_log10() +
  #scale_y_log10() +
  annotation_logticks(sides = "b")



raw_meta %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(nCount_RNA,percent.BCR+1)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  facet_wrap(.~ paste(Xp,sample)) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(sides = "b")
  
  
  
#~~~~~~~~~~~~~~~~~~~
# setting rough QC thresholds before centering across samples

ggplot() + annotate("text", x = 0, y = 0,
                    label = "setting rough QC thresholds before centering across samples") +
        theme_void()
        
raw_meta %>%
  ggplot(aes(nCount_RNA, paste(Xp,sample), fill = sample, color = Clog)) +
  geom_density_ridges(quantile_lines = TRUE) +
  scale_x_log10(limits = c(NA,NA)) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  geom_vline(xintercept = 1500)

raw_meta %>%
  ggplot(aes(nFeature_RNA, paste(Xp,sample), fill = sample, color = Clog)) +
  geom_density_ridges(quantile_lines = TRUE) +
  scale_x_log10(limits = c(NA,NA)) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  geom_vline(xintercept = 850)
  
  
  
  
#~~~~~~~~~~~~~~~~~~~
# centering and setting stringent cutoffs for low quality cells and doublet removal

lower_lim <- 1100
upper_lim <- 2600

test <- raw_meta %>%
  mutate(nFeat_raw = nFeature_RNA) %>%
  dplyr::filter(nFeature_RNA > 850 & !Clog & nCount_RNA > 1500) %>%
  mutate(nFeature_RNA = log10(nFeature_RNA)) %>%
  group_by(sample) %>%
  mutate(high = median(nFeature_RNA)) %>%
  ungroup() %>%
  mutate(nFeature_RNA = 10^(nFeature_RNA - high + mean(nFeature_RNA)))

test %>%
  ggplot(aes(nFeature_RNA, paste(Xp,sample), fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, bandwidth = 0.01) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_x_log10() +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  geom_vline(xintercept = c(lower_lim, upper_lim))



test %>%
  dplyr::filter(nFeature_RNA > lower_lim & nFeature_RNA < upper_lim & percent.mt < 2.5) %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(nCount_RNA,percent.mt)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  facet_wrap(.~ paste(Xp,sample)) +
  scale_x_log10() +
  #scale_y_log10() +
  annotation_logticks(sides = "b") +
  scale_y_continuous(limits = c(NA,10), oob = squish)


test %>%
  dplyr::filter(nFeature_RNA > lower_lim & nFeature_RNA < upper_lim & percent.mt < 2.5) %>%
  group_by(sample) %>%
  mutate(sample = paste(sample, "( n =", n(), ")")) %>%
ggplot(aes(nCount_RNA,nFeature_RNA)) +
  geom_hex(aes(fill = ..ndensity.., color = ..ndensity..),
           bins = 100) +
  geom_density2d(bins = 10, contour_var = "ndensity", color = "red") +
  facet_wrap(.~ paste(Xp,sample)) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(sides = "bl")


test %>%
  dplyr::filter(nFeature_RNA > lower_lim & nFeature_RNA < upper_lim & percent.mt < 2.5) %>%
  ggplot(aes(nFeature_RNA, paste(Xp,sample), fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, bandwidth = 0.005) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_x_log10() +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  geom_vline(xintercept = c(lower_lim, upper_lim))
  
  
#~~~~~~~~~~~~~~~~~
# amount of cells removed


kept <- test %>%
  dplyr::filter(nFeature_RNA > lower_lim & nFeature_RNA < upper_lim & percent.mt < 2.5) %>%
  pull(cell.id)



raw_meta <- raw_meta %>%
  mutate(keep = cell.id %in% kept)


ggplot(raw_meta, aes(sample, fill = keep, color = Clog)) +
  geom_bar() +
  theme_mrl() +
  facet_grid(.~Xp, scales = "free_x") +
  scale_fill_jco() +
  scale_color_aaas()


ggplot(raw_meta, aes(sample, fill = keep, color = Clog)) +
  geom_bar(position = "fill") +
  theme_mrl() +
  facet_grid(.~Xp, scales = "free_x") +
  scale_fill_jco() +
  scale_color_aaas()
  
  
  
  
dev.off()




#################################################
########### filter seurat object ################
#################################################
 
  

merged.seurat <- merged.seurat[,raw_meta$keep]


saveRDS(merged.seurat, file.path(results_dir,"seurat_objects","master_seurat_QCed.rds"))
  
  
  


