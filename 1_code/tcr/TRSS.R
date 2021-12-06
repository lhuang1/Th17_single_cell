##### Compare TCR similarities
## 1. EAE SPL_0/SPL_1 vs. other tissues
## 2. across all tissues
## 3. across all intra-tissue clusters
rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(dplyr)
library(tibble)
library(tidyr)

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()
  trt_gfp <- "EAE_GFPall" ## or "UT_GFPall"; treatment and GFP status combination
  clustering_date <- "2020-03-25"
  clonotyping_type_date <- "dominant_2020-04-02"
} else {
  today <- cargs[1]
  clustering_date <- cargs[2]
  clonotyping_date <- cargs[3]
}

dir_out <- paste0("2_pipeline/TCR/compare_tissue/")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

#### load data ####
tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL")
if (trt_gfp == "EAE_GFPall") {
  tissue_vec <- c(tissue_vec, "CNS", "DLN")
}
dir_cluster <- paste0("2_pipeline/clustering/", trt_gfp, "_intra/", clustering_date, "/")
dir_spl_cluster <- paste0("2_pipeline/clustering/", trt_gfp, "_SPL/", clustering_date, "/")

## meta data
meta <- list()
for (tissue in tissue_vec) {
  if (tissue == "SPL" & trt_gfp == "EAE_GFPall") {
    tmp  <- read.table(paste0(dir_spl_cluster, "FILES/meta_data_switched.csv"), header = T, row.names = 1, sep = ",")
    meta[[tissue]] <- tmp[,-grep("integrated_snn_res", colnames(tmp))] ## remove this column because the names are different across tissues
  } else {
    tmp  <- read.table(paste0(dir_cluster, "FILES/", tissue, "/meta_data.csv"), header = T, row.names = 1, sep = ",")
    meta[[tissue]] <- tmp[,-grep("integrated_snn_res", colnames(tmp))] ## remove this column because the names are different across tissues
  }
}
meta_all <- Reduce(rbind, meta)

## clonotypes
clty_all <- read.csv(paste0("2_pipeline/TCR/clonotype_assignment_", clonotyping_type_date, ".csv"), row.names = 1)

## merge clonotype info to meta data
meta_all$clonotype_id <- clty_all[rownames(meta_all),"clonotype_id"]
meta_all$pooled_clonotype_id <- paste(meta_all$treatment, meta_all$batch, meta_all$mouse, meta_all$clonotype_id, sep = "_")
meta_all$sample <- paste(meta_all$batch, meta_all$mouse, sep = "_") %>% factor()
meta_all$tissue_cluster <- factor(paste(meta_all$tissue, meta_all$seurat_clusters, sep = "_"))

meta_noNA <- meta_all[!is.na(meta_all$clonotype_id),] %>% droplevels()
stopifnot(all(meta_noNA$treatment == treatment)) ## check if treatment is correct
all_clonotypes <- unique(meta_noNA$pooled_clonotype_id)

if (trt_gfp == "EAE_GFPall") {
  ## similarity between SPL0/SPL1 to other tissues
  p_mtx_SPL01 <- sapply(c("SPL_0", "SPL_1"), function (x) {
    cluster_clonotypes <- meta_noNA$pooled_clonotype_id[meta_noNA$tissue_cluster == x]
    p <- table(cluster_clonotypes)[all_clonotypes]
    p[is.na(p)] <- 0
    p <- p / sum(p)
  })
  p_mtx_tissue <- sapply(tissue_vec, function (x) {
    cluster_clonotypes <- meta_noNA$pooled_clonotype_id[meta_noNA$tissue == x]
    p <- table(cluster_clonotypes)[all_clonotypes]
    p[is.na(p)] <- 0
    p <- p / sum(p)
  })
  sim_score <- sapply(c("SPL_0", "SPL_1"), function(x) {
    sapply(tissue_vec, function (y){
      sum(sqrt(p_mtx_SPL01[,x] * p_mtx_tissue[,y]))
    })
  })
  saveRDS(sim_score, file = paste0(dir_out, "similarity_score_SPL01_vs_tissue_", today, ".rds"))
}

## similarity between all pairs of tissues
sim_score <- sapply(tissue_vec, function(x) {
  sapply(tissue_vec, function (y){
    sum(sqrt(p_mtx_tissue[,x] * p_mtx_tissue[,y]))
  })
})
saveRDS(sim_score, file = paste0(dir_out, "similarity_score_tissue_", trt_gfp, "_", today, ".rds"))


## similarity between all pairs of tissue clusters
tissue_cluster_vec <- as.character(unique(meta_noNA$tissue_cluster))
p_mtx_tissue_cluster <- sapply(tissue_cluster_vec, function (x) {
  cluster_clonotypes <- meta_noNA$pooled_clonotype_id[meta_noNA$tissue_cluster == x]
  p <- table(cluster_clonotypes)[all_clonotypes]
  p[is.na(p)] <- 0
  p <- p / sum(p)
})
sim_score <- sapply(tissue_cluster_vec, function(x) {
  sapply(tissue_cluster_vec, function (y){
    sum(sqrt(p_mtx_tissue_cluster[,x] * p_mtx_tissue_cluster[,y]))
  })
})
saveRDS(sim_score, file = paste0(dir_out, "similarity_score_tissue_cluster_", trt_gfp, "_", today, ".rds"))

