#############################################################################################
########     UT GFPall Tregs (code modified based on make_figures/Fig2_efg.R)    ############
rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggsignif)
library(philentropy)
library(pheatmap)
library(grid)
set.seed(1)
source("1_code/utils.R")


#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-11-05"
  clustering_date <- "2020-03-25"
  assign_clonotype_date <- "2020-04-02"
} else {
  today <- cargs[1]
  clustering_date <- cargs[2]
  assign_clonotype_date <- cargs[3]
}

dir_out <- paste0("2_pipeline/TCR/Tregs/", today, "/")
if (!dir.exists(dir_out)) {dir.create(dir_out, recursive = T)}

### load data
cl_all <- read.csv(paste0("2_pipeline/TCR/clonotype_assignment_dominant_", assign_clonotype_date, ".csv"), row.names = 1, header = T)

tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL")
meta_ls <- sapply(tissue_vec, simplify = FALSE, FUN = function(tissue) {
  meta_tissue <- read.csv(paste0("2_pipeline/clustering/UT_GFPall_intra/", clustering_date, "/FILES/", tissue, "/meta_data.csv"), row.names = 1)
  ## merge clonotype with meta data (clustering)
  meta_tissue$clonotype_id <- cl_all[rownames(meta_tissue), "clonotype_id"]
  meta_tissue$pooled_clonotype_id <- paste(meta_tissue$orig.ident, meta_tissue$clonotype_id)
  
  #####################
  #*******debug@@@@@@@
  meta_tissue$batch <- as.character(meta_tissue$batch)
  meta_tissue$batch[meta_tissue$treatment == "UT" & meta_tissue$batch == "b8"] <- "b9"
  meta_tissue$batch[meta_tissue$treatment == "UT" & meta_tissue$batch == "b7"] <- "b8"
  meta_tissue$batch <- factor(meta_tissue$batch)
  
  meta_tissue
})



### percent treg/non-treg cells sharing TCR with SPL
### clonotype sharing
treg_cluster <- c(
  "MLN" = 3,
  "PP" = 3,
  "SI" = 2,
  "COL" = 2
)

mouse_vec <- c("m1", "m2")
n_shared_ls <- list()
for (tissue in names(treg_cluster)) {
  print(tissue)
  meta_ls[[tissue]]$is_treg <- meta_ls[[tissue]]$seurat_clusters == treg_cluster[tissue]
  meta_ls[[tissue]]$cluster_name <- paste0(tissue, "_", meta_ls[[tissue]]$seurat_clusters)
  if (tissue == "COL") {
    meta_ls[[tissue]]$cluster_name[meta_ls[[tissue]]$seurat_clusters %in% c(5, 6)] <- "COL_5&6" ## combine two proliferating clusters
  }
  cluster_name_sizes <- table(meta_ls[[tissue]]$cluster_name) %>% unlist()
  
  meta_tissue <- meta_ls[[tissue]]
  meta_spl <- meta_ls[["SPL"]]
  
  batch_vec <- meta_tissue$batch[meta_tissue$is_hashing] %>% droplevels() %>% levels()
  n_shared <- data.frame(matrix(NA, nrow = 4, ncol = 6))
  colnames(n_shared) <- c("batch", "mouse", "Treg_total", "Treg_shared", "nTreg_total", "nTreg_shared")
  meta_tissue_by_sample <- list()
  i = 1
  for (batch in batch_vec) {
    for (mouse in mouse_vec) {
      sample_name <- paste(batch, mouse, sep = "_")
      cat(sample_name, "\n")
      meta_tissue_sub <- meta_tissue[meta_tissue$batch == batch & meta_tissue$mouse == mouse & !is.na(meta_tissue$clonotype_id),] %>% droplevels()
      meta_spl_sub <- meta_spl[meta_spl$batch == batch & meta_spl$mouse == mouse & !is.na(meta_spl$clonotype_id),] %>% droplevels()
      meta_tissue_by_sample[[sample_name]] <- meta_tissue_sub
      n_treg <- sum(meta_tissue_sub$is_treg)
      n_ntreg <- sum(!meta_tissue_sub$is_treg)
      n_treg_shared <- sum(meta_tissue_sub$clonotype_id[meta_tissue_sub$is_treg] %in% meta_spl_sub$clonotype_id)
      n_ntreg_shared <- sum(meta_tissue_sub$clonotype_id[!meta_tissue_sub$is_treg] %in% meta_spl_sub$clonotype_id)
      n_shared[i,"batch"] <- batch
      n_shared[i,"mouse"] <- mouse
      n_shared[i,"Treg_total"] <- n_treg
      n_shared[i,"Treg_shared"] <- n_treg_shared
      n_shared[i,"nTreg_total"] <- n_ntreg
      n_shared[i,"nTreg_shared"] <- n_ntreg_shared
      i <- i + 1
    }
  }
  n_shared$Treg_frac <- n_shared$Treg_shared / n_shared$Treg_total
  n_shared$nTreg_frac <- n_shared$nTreg_shared / n_shared$nTreg_total
  
  for (i in seq_along(n_shared[,1])) {
    print(i)
    if (all(is.na(n_shared[i, c("Treg_total", "nTreg_total")])) | all(n_shared[i, c("Treg_total", "nTreg_total")] == 0)) {
      n_shared[i,"p_chisq"] <- NA
      n_shared[i,"p_fisher"] <- NA
    } else {
      test_chisq <- chisq.test(x = matrix(data = as.numeric(n_shared[i, 3:6]), nrow = 2))
      test_fisher <- fisher.test(x = matrix(data = as.numeric(n_shared[i, 3:6]), nrow = 2))
      n_shared[i,"p_chisq"] <- test_chisq$p.value
      n_shared[i,"p_fisher"] <- test_fisher$p.value
    }
  }
  n_shared_ls[[tissue]] <- n_shared
}

saveRDS(n_shared_ls, file = paste0(dir_out, "treg_share_with_SPL_", today, ".rds"))


