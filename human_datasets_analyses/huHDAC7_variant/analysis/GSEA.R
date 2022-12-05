#####################################################################################################################
############################################## Set up ###############################################################
#####################################################################################################################


# Import packages and custom wrapper functions
source("src/bulk_wrappers.R")

# other libraries
library(fgsea)
library(AnnotationDbi)
library(eulerr)
library(DOSE)
library(proxy)
library(reshape2)
library(ComplexHeatmap)
library(ggdendro)
library(pheatmap)
library(gtable)



#####################################################################################################################
########################################### Prepare preranked list ##################################################
#####################################################################################################################


#Import DESeq2 output
degDESeqWT <- read.delim("huHDAC7_variant/results/DGE/degDESeqWT.txt")
degDESeqMut2 <- read.delim("huHDAC7_variant/results/DGE/degDESeqMut.txt")
degDESeqVariant <- read.delim("huHDAC7_variant/results/DGE/degDESeqVariant.txt")

# GSEA requirements

# Based on:
#   
#   https://stephenturner.github.io/deseq-to-fgsea/
#   
#   http://evomics.org/wp-content/uploads/2018/10/RNASeqWorkFlow-1.html


#Remove potential duplicates entries based on gene names
#(they can be unprocessed transcript or alternative mapping to neighbor locus)
#We arbitrarily choose the entry that has the highest baseMean expression.
#It is advised to use gene symbols for input in GSEA
#I checked the difference between gene symbol and gene name and there is no difference.
#So using the gene_name column should be ok.
degDESeqVariant.nodup <- degDESeqVariant %>% 
  group_by(gene_name) %>% 
  top_n(n = 1, wt = baseMean) %>%
  ungroup()

degDESeqWT.nodup <- degDESeqWT %>% 
  group_by(gene_name) %>% 
  top_n(n = 1, wt = baseMean) %>%
  ungroup()

degDESeqMut2.nodup <- degDESeqMut2 %>% 
  group_by(gene_name) %>% 
  top_n(n = 1, wt = baseMean) %>%
  ungroup()


#Generate the ranked file by selecting the gene name and the statistics variables and remove NAs
RNKvariant <- degDESeqVariant.nodup %>%
  mutate(Gene_Name = gene_name) %>%
  dplyr::select(Gene_Name, stat) %>%
  dplyr::filter(is.na(Gene_Name)==FALSE)




#####################################################################################################################
########################################### GSEA on hallmarks #######################################################
#####################################################################################################################

#########################
## Analysis with fgsea ##
#########################
#Load hallmark GSEA genesets:
pathways.hallmark <- gmtPathways("huHDAC7_variant/data/h.all.v6.2.symbols.gmt")

#Convert variant results to a vector:
RNKvariantVector <- deframe(RNKvariant)

#Run fgsea on variant set
fgseaVarH <- fgsea(pathways=pathways.hallmark, stats=RNKvariantVector, nperm=50000, minSize = 50)

topPathwaysUp <- fgseaVarH[ES > 0][head(order(-NES), n=10), pathway]
topPathwaysDown <- fgseaVarH[ES < 0][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.hallmark[topPathways], RNKvariantVector, fgseaVarH, 
              gseaParam = 0.5, c(5,5,1,1,1))


########################
## Analysis with DOSE ##
########################

# Same analysis with the DOSE wrapper that provides a nice interface for making enrichment plots
# the gseaplot2 from enrichplot package not available on R3.5 so I sourced the function.

detach("package:dplyr", unload = TRUE)
library(clusterProfiler)
library(dplyr)
library(enrichplot)
source("src/gseaplot2_PPAmod.R")

RNK_H7 <- degDESeqVariant.nodup %>%
  mutate(Gene_Name = gene_name) %>%
  dplyr::select(Gene_Name, stat) %>%
  dplyr::filter(is.na(Gene_Name) == FALSE) %>%
  arrange(-stat)

#Convert variant results to a vector:
RNK_H7Vector <- deframe(RNK_H7)

hallmark.df <- list_to_df(pathways.hallmark) %>%
  unnest() %>%
  plyr::rename(c("name" = "Term")) %>%
  plyr::rename(c("list.element" = "Gene"))


gse <-
  GSEA(
    geneList = RNK_H7Vector,
    nPerm = 50000,
    minGSSize = 50,
    maxGSSize = 5000,
    TERM2GENE = hallmark.df,
    pvalueCutoff = 0.1
  )


# plot with OXPHOS enrichment

gseaplot2(
  gse,
  geneSetID = 1,
  title = gse$Description[2],
  color = "black"
)


ggsave("huHDAC7_variant/results/figs/gseaplots.pdf",
       width = 6,
       height = 5.5)



# plot with significant genesets


gseResult <- gse@result %>%
  separate(leading_edge,
           into = c("tags", "list", "signal"),
           sep = ",") %>%
  mutate(Overlap = as.integer(gsub("[^0-9.-]", "", tags))) %>%
  mutate(list = as.integer(gsub("[^0-9.-]", "", list))) %>%
  mutate(signal = as.integer(gsub("[^0-9.-]", "", signal)))


ggplot(gseResult %>% dplyr::mutate(ID = gsub("_", " ", ID)),
       aes(
         x = ID,
         y = NES,
         size = Overlap,
         fill = -log10(pvalue)
       )) +
  geom_point(shape = 21,
             stroke = 0.1,
             color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-max(abs(gseResult$NES)), max(abs(gseResult$NES)))) +
  coord_flip() +
  scale_size_continuous(range = c(1, 6)) +
  scale_fill_viridis(limits = c(1.30103, NA)) +
  labs(x = "", y = "Normalized enrichment score") +
  scale_x_discrete(
    position = "bottom",
    labels = function(x)
      str_wrap(x, width = 12)
  ) +
  theme_mrl(0.6) +
  theme(
    axis.title.y = element_blank(),
    axis.ticks = element_line(),
    axis.line = element_blank(),
    axis.text = element_text(colour = "black"),
    legend.position = "right",
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    aspect.ratio = 1,
    legend.key.size = unit(0.3, "cm")
  )





ggsave("huHDAC7_variant/results/figs/gseaplots2.pdf",
       width = 5,
       height = 4.5)




#####################################################################################################################
########################################### GSEA on HDAC7 KI ########################################################
#####################################################################################################################


clustered_DEGs <- readRDS("../msHDAC7KI_analyses/tables/Treg_poisson_glm_DGE_results.csv") %>%
  dplyr::filter(fdr < 0.05 & !grepl("^TRB|^TRA|^RPS|^RPL|^MT-", gene_name, ignore.case = T))


correspondance <- convertMouseGeneList(clustered_DEGs$gene_id) %>%
  dplyr::select(Gene.stable.ID, Gene.name) %>% 
  dplyr::filter(!is.na(Gene.name)) %>%
  distinct()

clustered_DEGs <- clustered_DEGs %>%
  left_join(correspondance, by = c("gene_id" = "Gene.stable.ID")) %>%
  group_by(Gene.name) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  dplyr::filter(count == 1 | Gene.name == toupper(external_gene_name))


# gsea on clusters
clustered_DEGs_list <- clustered_DEGs %>% 
  dplyr::select(cluster, Gene.name) %>%
  mutate(cluster = paste0("cluster",cluster)) %>%
  distinct() %>%
  group_by(cluster) %>%
  nest() %>%
  mutate(data = purrr::map(data, ~pull(.,Gene.name))) %>%
  deframe()



msKI <- list_to_df(clustered_DEGs_list) %>% 
  unnest() %>% 
  plyr::rename(c("name" = "Term")) %>% 
  plyr::rename(c("list.element" = "Gene"))


gse <- GSEA(
  geneList = sort(RNKvariantVector, decreasing = T),
  nPerm = 50000,
  minGSSize = 50,
  maxGSSize = 5000,
  TERM2GENE = msKI,
  pvalueCutoff = 0.1
)


title <- gse@result %>%
  dplyr::filter(ID == "upKI") %>%
  mutate(NES = round(NES, 2)) %>%
  #mutate(pvalue = as.character(format(pvalue))) %>%
  dplyr::select(NES, pvalue) %>%
  gather() %>%
  mutate(value = signif(value, 3)) %>%
  unite(merged, c(key, value), sep = " = ") %>%
  pull(merged) %>%
  paste(collapse = ", ") %>%
  paste("HDAC7 KI brain infiltrating Treg geneset.", .)


gseaplot2(
  gse,
  geneSetID = "upKI",
  title = title,
  color = "black",
  base_size = 3
) 

ggsave("huHDAC7_variant/results/figs/huVar_vs_msVar.pdf", width=3, height = 2)


