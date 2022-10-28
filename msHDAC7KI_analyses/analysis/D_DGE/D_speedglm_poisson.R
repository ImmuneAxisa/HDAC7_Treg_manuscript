
## For R3.5!

# This is from v6 of speedglm benchmarking
# modification: run parallel iterations without creating a dense matrix
# instead transpose the GEX matrix as subsetting columns (genes) is faster: about 1sec
# Still should create an extra ~4hours at 13,000 genes * 1sec
# but might be better for memory usage

# Extra modification from v5: dumb here to have more than one chunk per worker
# so create n subset of genes to parallelise that equal n number of workers

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
library(broom)
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
library(future.apply)
library(speedglm)
library(rslurm)



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

# function for building chunk of genes for //
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


#########################################################
################### variables ###########################
#########################################################

# parallelization parameters:
nodes = 1
cpus_per_node = 8





#########################################################
################### paths and stems #####################
#########################################################


# Set wd on project folder
setwd("/home/pa326/project/HDAC7_KI_EAE_10x/")


results_dir <- "results/master"
dir.create(results_dir, showWarnings = FALSE)


objects_dir <- file.path(results_dir,"seurat_objects")
dir.create(objects_dir, showWarnings = FALSE)



results_dir <- file.path(results_dir,"DGE")
dir.create(results_dir, showWarnings = FALSE)

#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- readRDS(file.path(objects_dir,"master_seurat_harmony_Tcells_low_res_clustered.rds"))


#########################################################
################ Select genes to test ###################
#########################################################





genes_to_test <- merged.seurat@assays[["RNA"]]@meta.features %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("gene_name") %>%
                    dplyr::filter(vst.mean > 10^-3) %>%
                    pull(gene_name)
                    
                    
#genes_to_test <- genes_to_test[1:1000]


cat("\nnumber of genes to test \n",length(genes_to_test))

GEX <- merged.seurat@assays$RNA@counts

# binarize
# and transpose because it's faster to subset columns
GEX <- t(GEX)

#########################################################
############### Get MELD probabilities ##################
#########################################################


meta <- merged.seurat@meta.data %>%
  tibble::rownames_to_column("index") %>%
  mutate(SNNk20_algo1_res0.3 = as.character(SNNk20_algo1_res0.3)) %>%
  mutate(Sub_res0.15 = as.character(Sub_res0.15)) %>%
  mutate(Sub_res0.1 = as.character(Sub_res0.1)) %>%
  mutate(curated_cluster = ifelse(SNNk20_algo1_res0.3 %in% c("0","1"), Sub_res0.15, Sub_res0.1)) %>%
  mutate(it_all = TRUE) 


    
        
cat("\nnumber of cells initially:", 
    nrow(meta),
    "\n\n")
            
if(identical(rownames(GEX), meld$index))   {
    cat("\n cell orderings are concordant \n")
}  else {
    cat("\n GEX matrix and MELD prob cell orderings do not match \n")
    quit(save = "no", status = 1)
}

rm(merged.seurat)

############################################################################ 
############################################################################
############## Prepare DGE analysis ########################################
############################################################################
############################################################################ 



###########################################################################################
##################### Parallelise chunks for DGE testing ##################################
########################################################################################### 


cat("\n start job submission \n")

# define param for export of results

output_prefix <- file.path(results_dir,"poisson_glm_labelPred_")
pdf(paste0(output_prefix, "cpc.pdf"))

# this will set all temp files for slurm submissions to be written to scratch
setwd("/gpfs/ycga/scratch60/hafler/pa326/HDAC7_KI_EAE_10x/")


# define which iteration and cluster to run DGE on
groupings <- grep("^it_", colnames(meld), value = T)
clusters_to_test <- unique(meld$curated_cluster)
clusters_to_test <- grep("0_|4|6", clusters_to_test, value = T)


#--------------------------#--------------------------#--------------------------#--------------------------
# iteration for each predictor
for(group in groupings){

    # term for formula
    #pred_term <- pred
    pred_term <- "Geno"
    
    form <- paste("genes ~", pred_term, " + nCount_RNA")
    
    cat("\n",form, "\n")
    
    #--------------------------#--------------------------#--------------------------#--------------------------
    # Iteration on each cluster
    for(clust in clusters_to_test) {
        
        cat("\nrunning", clust, "on", group, "combination\n")
        #--------------------------#--------------------------#--------------------------
        # subset the data to cluster and appropriate patient combination
        tmp_meld <- meld %>%
            dplyr::filter(!!sym(group)) %>%
            dplyr::filter(curated_cluster == clust)
            
            
            
        cat("number of cell:", nrow(tmp_meld), "\n")
        
        tmp_GEX <- GEX[tmp_meld$index,]
        
        #--------------------------#--------------------------#--------------------------
        # filter out low expression genes.
        
        cpc <- Matrix::colMeans(tmp_GEX)
        
        
        plot <- qplot(x = cpc + 10^-4) + 
            geom_vline(xintercept = 0.005, color = "red") + 
            labs(title = paste("poisson_glm on cluster", clust, "with", nrow(tmp_GEX), "cells \n", group, ":", form),
                    caption = paste(sum(cpc > 0.005), "genes pass of", length(cpc))) +
            scale_x_log10()
        
        
        
        print(plot_grid(plot))
            
        cat("\nnumber of genes before filtering:",
            ncol(tmp_GEX), 
            "\n")
        
        tmp_GEX <- tmp_GEX[,cpc > 0.005]
        
        cat("\nnumber of genes after filtering:",
            ncol(tmp_GEX), 
            "\n")
            
            
        chunks <- chunk2(colnames(tmp_GEX), nodes * cpus_per_node)
        
        
        cat(length(chunks), "gene chunks to split on", nodes, "nodes with", cpus_per_node, "cpus_per_node", 
            sep = "\n")
            
        

        # prepare slurm jobs with rslurm:
        
        
        slr_job <- slurm_map(chunks, function(x) {
        
        
            # start timer for GEX matrix subsetting and bringing to dense matrix
            start_time <- Sys.time()
            
            # subset main GEX matrix for given chunk of genes
            #GEX2 <- as.matrix(GEX[,x])
            
            message(head(x))
            print((head(x)))
            
            
            
            if(T) {
            
            #end_time <- Sys.time()
            #message("time to prepare subset the GEX for 1 chunk")
            #message(capture.output(end_time - start_time))
            
            
            #message(paste("start chunk with", length(x), "genes"))
            
            #x <- x[1:70] # only evaluate 10% of genes for benchmarking
            
            # Run correlation on the chunk
            results <- list()
            for(genes in x) {
                
                
                # start a timer
                start_time <- Sys.time()
                
                print(genes)
                print(start_time)
                print(tmp_meld[1:2,1:4])
                
                
                # prepare data frame for regression
                logit.frame <- data.frame(genes = tmp_GEX[,genes]) %>%
                    tibble::rownames_to_column("index") %>%
                    left_join(tmp_meld, by = "index")
                    
                    
                # run logit with speedglm
                model <- tryCatch(speedglm(as.formula(form), 
                                        data = logit.frame, 
                                        family = poisson()),
                                error = capture.output)
                                
                print(difftime(Sys.time(), start_time, units = "secs"))
                
                # add useful metadata for later
                model[["gene_name"]] <- genes
                model[["runtime"]] <- difftime(Sys.time(), start_time, units = "secs")
                                        
                
                
                # add to the results object
                results[[genes]] <- model
                
                
                # every hundred genes print progress
                if(is.integer(length(results) / 100)) {
                    message(paste(length(results), "genes processed\n"))
                }
                rm(logit.frame, model)
                gc()
                }
            
            gc()
            return(results)
            
            }
            
            
            },
            # slurm options
            jobname = paste(output_prefix,group,clust,
                        gsub(" |:", "_",date()),
                        sep = "_"),  #add date to append to tmp directories created (it avoids rewriting and mixing files with previous iterations)
            nodes = nodes, 
            cpus_per_node = cpus_per_node,
            preschedule_cores = FALSE,
            global_objects = c("tmp_GEX","tmp_meld", "form"),
            slurm_options = list(mem = "32gb", 
                                partition = "general", 
                                error = "slurm_%a.err",
                                time = "12:00:00"),
            submit = T
            )
            
    ############################################################################ 
    ############################################################################ 
    ##################### Retrieve results and save output ##################### 
    ############################################################################ 
    ############################################################################ 
        
        
        
        
        
        print("finished array submission, submit data retrieval job with dependency")
        
        slurm_call(function() {
        
        
                setwd("/gpfs/ycga/scratch60/hafler/pa326/HDAC7_KI_EAE_10x/")
        
            
                # get results and logs
                ######################
                    
                logs <- get_job_status(slr_job)
                
                res <- get_slurm_out(slr_job)
                
                
                # write results to file
                #######################
                
                setwd("/home/pa326/project/HDAC7_KI_EAE_10x/")
                
                #results_dir <- "results/DGE/"
                
                saveRDS(logs, paste0(output_prefix,paste(group,clust,"logs.rds",sep = "_")))
                
                saveRDS(res, paste0(output_prefix,paste(group,clust,"res.rds",sep = "_")))
                
                
                # summarize lm object data
                ##########################
                
                results <- do.call(c,res)
                
                rm(res)
        
                print("list length after collapsing embedded lists")
                print(length(results))
        
                if(T) {
                
                metadata <- lapply(results, function(model) {
                    
                    if(length(model) < 4) {
                        model[["convergence"]] <- F
                        model[["aic"]] <- model[[1]] # put error message in place of AIC if failed
                    }
                    
                
                    data.frame(
                                gene_name = model[["gene_name"]],
                                runtime = as.double(model[["runtime"]]), 
                                converged = model[["convergence"]],  ## this slot is called convergence for speedglm and converged for base R glm. Confusing!
                                aic = model[["aic"]]
                            )
                
                })
                
                metadata <- do.call(rbind,metadata)
                
                metadata <- metadata %>%
                    dplyr::left_join(data.frame("cpc" = cpc) %>% tibble::rownames_to_column("gene_name"))
                
                saveRDS(metadata, paste0(output_prefix,paste(group,clust,"meta.rds",sep = "_")))
                
                }
                
                #########################################################
                ##################### get summary #######################
                #########################################################
                
                summarized_results <- lapply(results, function(model) {
                    
                    if(length(model) > 3) {
                    
                    model %>%
                        broom::tidy() %>%
                        mutate(gene_name = model[["gene_name"]])
                
                    }  else {
                        return(NULL)
                    } 
                                                
                })
        
                summarized_results <- do.call(rbind,summarized_results)
            
                saveRDS(summarized_results, paste0(output_prefix,paste(group,clust,"summary.rds",sep = "_")))
            
                
                
            
            
                },
            global_objects = c("slr_job", "group","clust", "output_prefix", "cpc"),
            jobname = paste("speedglm_gather",group,clust, sep = "_"),
            slurm_options = list(mem = "16gb", 
                                partition = "general", 
                                error = "slurm_%a.err",
                                dependency=sprintf("afterany:%s", slr_job$jobid),
                                time = "02:00:00"),
            submit = T
            )
    

    }
}

dev.off()
