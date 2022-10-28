#####################################################################################################################
############################################## Set up ###############################################################
#####################################################################################################################


##############
## packages ##
##############

library(knitr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(viridis)
library(scales)
library(ggrepel)
library(tximport)
library(edgeR)
library(limma)
library(DESeq2)
library(RUVSeq)
library(circlize)
library(Hmisc)
library(DEFormats)
library(statmod)
library(vsn)
library(gplots)
library(cowplot)
library(ggsci)
library(WGCNA)
library(cluster)
library(outliers)
library(gtable)
library(data.table)
library(ggExtra)
library(gridGraphics)

options(stringsAsFactors = FALSE)




##################
## ggplot theme ##
##################

theme_mrl <- function(x = 1) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(size = 12*x,face = "bold", angle = 0, vjust = 0.6),
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

#####################################################################################################################
#################################### Import quantitation results ####################################################
#####################################################################################################################


# Data from GSE209596
treg.counts <- read.delim("huHDAC7_coexpression/results/treg.rsem.filtered.counts.txt") %>%
  dplyr::mutate(gene_anno = paste0(gene_id,"_",gene_name)) %>%
  tibble::column_to_rownames("gene_anno") %>%
  dplyr::select(-gene_id, -gene_name) %>%
  as.matrix()
treg.counts <- cpm(treg.counts, log = F)


treg.counts.Exvivo <- treg.counts[,grepl("Exvivo",colnames(treg.counts))] 
treg.counts.Exvivo <- treg.counts.Exvivo[,!grepl("CSF",colnames(treg.counts.Exvivo))] 




#####################################################################################################################
##################################### Build coexpression network ####################################################
#####################################################################################################################

# From https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-05-NetworkConstruction.pdf

datExpr=t(cpm(treg.counts.Exvivo, log = T))

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red", lty = 3)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")




## here we define the adjacency matrix using soft thresholding with beta=6
#ADJ1=abs(cor(t(treg.counts.Exvivo),use="p"))^6
## When you have relatively few genes (<5000) use the following code
# k=as.vector(apply(ADJ1,2,sum, na.rm=T))
## When you have a lot of genes use the following code
k=softConnectivity(datE=t(cpm(treg.counts.Exvivo, log = T)),power=6)
## Plot a histogram of k and a scale free topology plot
#sizeGrWindow(10,5)
#par(mfrow=c(1,2))
#hist(k)
#scaleFreePlot(k, main="Check scale free topology\n")


#restrict:
datExpr=t(cpm(treg.counts.Exvivo, log = T))[, rank(-k,ties.method="first" )<=10000]

# Choose parameters for deepSplit and mergeCutHeight
x <- c(0)
y <- c(0.11)

for (deepSplit in x) {
  for(mergeCutHeight in y) {
    ##############################
    # one step module construction
    
    net = blockwiseModules(datExpr, power = 6,
                           TOMType = "signed Nowick", minModuleSize = 30,
                           reassignThreshold = 1e-1000, mergeCutHeight = mergeCutHeight,
                           deepSplit = deepSplit,
                           detectCutHeight = 0.95,
                           numericLabels = TRUE, pamRespectsDendro = T,
                           saveTOMs = F,
                           maxBlockSize = 10000,
                           verbose = 3)
    
    
    
    # open a graphics window
    sizeGrWindow(12, 9)
    # Convert labels to colors for plotting
    Colors = labels2colors(net$unmergedColors)[net$blockGenes[[1]]]
    mergedColors = labels2colors(net$colors)[net$blockGenes[[1]]]
    
    test <- data.frame(mod = net$colors) %>%
      tibble::rownames_to_column("genes")
    
    HDAC7module <- test %>% dplyr::filter(grepl("HDAC7", genes)) %>% dplyr::pull(mod)
    test <- test %>% dplyr::filter(mod == HDAC7module)
    
    HDAC7colors = labels2colors(net$colors == HDAC7module)[net$blockGenes[[1]]]
    
    # Plot the dendrogram and the module colors underneath
    #pdf(file="../results/MattRNAseqCounts/WGCNAdendrogramTregs.pdf")
    plotDendroAndColors(net$dendrograms[[1]], 
                        cbind(Colors, mergedColors, HDAC7colors),
                        c("Unmerged", "Merged", "HDAC7 module"),
                        dendroLabels = FALSE, hang = 0.03, 
                        addGuide = TRUE, guideHang = 0.05,
                        main = paste0("deepSplit=",
                                      deepSplit,
                                      ", ",
                                      "mergeCutHeight=", 
                                      mergeCutHeight))
    #dev.off()
    
    
    collectGarbage()
    
    
  }
}
#####################################################################################################################
######################################## save results to file #######################################################
#####################################################################################################################



#  tweak colors and save and save dendrogram
sizeGrWindow(2, 9)
mergedColors = labels2colors(net$colors, 
                             colorSeq = colorRampPalette(palette(pal_jco(alpha=1)(10)))(24))[net$blockGenes[[1]]]
pdf(file="huHDAC7_coexpression/results/figs/WGCNAdendrogramTregsLarge.pdf", width = 10, height = 7)
plotDendroAndColors(net$dendrograms[[1]], 
                    cbind(mergedColors),
                    c("Modules"),
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05,
                    main = "")
dev.off()
length(unique(mergedColors))



#plot correspondance to module numbers
mergedColors = labels2colors(net$colors, 
                             colorSeq = colorRampPalette(palette(pal_jco(alpha=1)(10)))(24))[net$blockGenes[[1]]]
module.names <- c("palette")
for(i in unique(net$colors)) {
  mergedColors = cbind(mergedColors,labels2colors(net$colors == i))
  module.names <- c(module.names, i) 
}
plotDendroAndColors(net$dendrograms[[1]], 
                    cbind(mergedColors),
                    module.names,
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05,
                    main = "")


# document color correspondance to module number 
module.colors <- data.frame(Module = net$colors, color = mergedColors) %>%
  unique() %>%
  arrange(Module) %>%
  dplyr::filter(Module != 0) %>%
  mutate(Module = as.factor(Module))

ggplot(module.colors, aes(x = 1, y = Module, fill = Module)) +
  geom_tile() +
  scale_fill_manual(values = module.colors$color) +
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        aspect.ratio = 20)

ggsave("huHDAC7_coexpression/results/figs/Module_colors.pdf", width = 3, height = 20)

# Number of genes per module
module.counts <- data.frame(Module = net$colors) %>%
  group_by(Module) %>%
  summarise(genes = n())

# Export HDAC7 module genes
H7mod <- test %>%
  separate(genes, c("gene_id", "gene_name"), sep = "_") %>%
  dplyr::pull(gene_name)
write.table(test %>% separate(genes, c("gene_id", "gene_name"), sep = "_"),
            "huHDAC7_coexpression/results/MattRNAseqCounts/HDAC7moduleTreg_wID.txt", row.names = F, col.names = T, quote = F)
write.table(H7mod, "huHDAC7_coexpression/results/MattRNAseqCounts/HDAC7moduleTreg.txt", row.names = F, col.names = F, quote = F)


#####################################################################################################################
####################################### Enrichr pathways plots ######################################################
#####################################################################################################################

# enrichr analyses were performed on the enrichr webportal using the HDAC7 module as input

###############
### TF PPIs ###
###############

PPIs <- read.delim("huHDAC7_coexpression/results/enrichr/Transcription_Factor_PPIs_table.txt",stringsAsFactors = F) %>%
  mutate(Rep = str_count(Genes, ";")+1) %>%
  top_n(10, Combined.Score) %>%
  mutate(Term = as.factor(Term)) %>%
  mutate(Term = forcats::fct_reorder(Term, Combined.Score)) %>%
  separate(Overlap, sep = "/", c("N","D"), convert = T) %>%
  mutate(Overlap = round(N/D*100*log2(N+1))) %>%
  dplyr::select(-N,-D)



ggplot(PPIs, aes(x = Term , y = Combined.Score, fill = -log10(Adjusted.P.value), size = Rep)) +
  geom_point(shape = 21, stroke = 0.1) +
  coord_flip() +
  theme_mrl(0.8) +
  theme(axis.title.y=element_blank(),
        axis.ticks = element_line(),
        axis.line = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        aspect.ratio = 1,
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(position = "bottom",labels = function(x) str_wrap(x, width = 50)) +
  labs(y = "Combined Score") +
  scale_fill_viridis() +
  guides(size= guide_legend(title = "Number of interactants"))

ggsave("huHDAC7_coexpression/results/enrichr/PPIenrichment.pdf", width = 5, height = 3)



#############
### GO CC ###
#############


GO_CC <- read.delim("huHDAC7_coexpression/results/enrichr/GO_Cellular_Component_2018_table.txt",stringsAsFactors = F) %>%
  mutate(Rep = str_count(Genes, ";")+1) %>%
  top_n(10, Combined.Score) %>%
  dplyr::filter(Adjusted.P.value < 0.1) %>%
  separate(Term, sep = "\\(GO", c("Term"), drop = "drop") %>%
  mutate(Term = as.factor(Term)) %>%
  mutate(Term = forcats::fct_reorder(Term, Combined.Score)) %>%
  separate(Overlap, sep = "/", c("N","D")) %>%
  mutate(Overlap = round(as.integer(N)/as.integer(D)*100)) %>%
  dplyr::select(-N,-D)


ggplot(GO_CC, aes(x = Term , y = Combined.Score, fill = -log10(Adjusted.P.value), size = Rep)) +
  geom_point(shape = 21, stroke = 0.1) +
  coord_flip() +
  theme_mrl(0.6) +
  theme(axis.title.y=element_blank(),
        axis.ticks = element_line(),
        axis.line = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        aspect.ratio = 1,
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(position = "bottom",labels = function(x) str_wrap(x, width = 30)) +
  labs(y = "Combined Score") +
  scale_fill_viridis() +
  guides(size= guide_legend(title = "Number of genes"))

ggsave("huHDAC7_coexpression/results/figs/GO_CC.pdf", width = 5, height = 3)
