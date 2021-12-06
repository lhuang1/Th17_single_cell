#######################################################
######## TCR sharing between Tregs and SPL ############
#######################################################
## Compute proportion of treg/non-treg cells in each tissue sharing TCR with SPL (UT, GFP all)

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(dplyr)
library(tidyr)
library(openxlsx)


#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()
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
  meta_tissue$pooled_clonotype_id <- paste(meta_tissue$batch, meta_tissue$mouse, meta_tissue$clonotype_id, sep = "_")
  meta_tissue[!is.na(meta_tissue$raw_clonotype_id),]
})

## use clsuter annotation to filter out non-conventional Th17 cells (Treg-like and proliferating clusters)
cluster_ann <- read.xlsx(paste0("2_pipeline/clustering/UT_GFPall_intra/", clustering_date, "/cluster_annotation.xlsx"), sheet = "processed")
cluster_ann$tissue <- sub("_.*", "", cluster_ann$Cluster)
tmp <- cluster_ann$Cluster[which(cluster_ann$Annotation == "Treg-like")]
treg_cluster <- as.integer(sub(".*_", "", tmp)) %>% `names<-`(sub("_.*", "", tmp))



## compare TCR
SPL_clonotypes <- unique(meta_ls$SPL$pooled_clonotype_id)
n_shared_ls <- list()
for (tissue in names(treg_cluster)) {
  print(tissue)
  meta_ls[[tissue]]$cluster_name <- paste0(tissue, "_", meta_ls[[tissue]]$seurat_clusters)
  if (tissue == "COL") {
    meta_ls[[tissue]]$cluster_name[meta_ls[[tissue]]$seurat_clusters %in% c(5, 6)] <- "COL_5&6" ## combine two proliferating clusters
  }
  # cluster_name_sizes <- table(meta_ls[[tissue]]$cluster_name) %>% unlist()
  # colnames(n_shared) <- c("batch", "mouse", "Treg_total", "Treg_shared", "nTreg_total", "nTreg_shared")
  n_shared <- meta_ls[[tissue]] %>% 
    mutate(is_treg = seurat_clusters == treg_cluster[tissue]) %>% 
    group_by(batch, mouse, .drop = FALSE) %>% 
    mutate(shared_SPL = pooled_clonotype_id %in% SPL_clonotypes) %>% 
    summarise(Treg_total = sum(is_treg),
              Treg_shared = sum(is_treg & shared_SPL),
              nTreg_total = sum(!is_treg),
              nTreg_shared = sum((!is_treg) & shared_SPL),
              Treg_frac = Treg_shared / Treg_total,
              nTreg_frac = nTreg_shared / nTreg_total) %>% 
    data.frame()
  ## test for equal proportions
  for (i in seq_along(n_shared[,1])) {
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


