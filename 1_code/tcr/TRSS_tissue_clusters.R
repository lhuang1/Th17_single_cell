rm(list = setdiff(ls(), "so_all"))

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(gridExtra)
set.seed(1)
source("1_code/utils.R")
source("1_code/tcr/utils.R")

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-11-02"
  clustering_date <- "2020-03-25"
  clonotyping_type_date <- "dominant_2020-04-02"
  prep_type_date <- "dominant_TCR_2020-03-24"
} else {
  today <- cargs[1]
  clustering_date <- cargs[2]
  clonotyping_date <- cargs[3]
}

dir_out <- paste0("2_pipeline/TCR/compare_tissue/")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

#### load data ####
treatment <- "EAE"
tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL", "CNS", "DLN")
dir_cluster <- paste0("2_pipeline/clustering/EAE_GFPall_intra/", clustering_date, "/")
dir_spl_cluster <- paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/")

## meta data
meta <- list()
for (tissue in tissue_vec) {
  if (tissue == "SPL" & treatment == "EAE") {
    tmp  <- read.table(paste0(dir_spl_cluster, "FILES/meta_data_switched.txt"), header = T, sep = " ")
    meta[[tissue]] <- tmp[,-grep("integrated_snn_res", colnames(tmp))]
  } else {
    tmp  <- read.table(paste0(dir_cluster, "FILES/", tissue, "/meta_data.csv"), header = T, row.names = 1, sep = ",")
    meta[[tissue]] <- tmp[,-grep("integrated_snn_res", colnames(tmp))]
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
stopifnot(all(meta_noNA$treatment == treatment)) ## check if treatment is correct (should be EAE cells only)

## similarity to tissues
all_clonotypes <- unique(meta_noNA$pooled_clonotype_id)
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

## similarity between all pairs of tissues
sim_score <- sapply(tissue_vec, function(x) {
  sapply(tissue_vec, function (y){
    sum(sqrt(p_mtx_tissue[,x] * p_mtx_tissue[,y]))
  })
})
saveRDS(sim_score, file = paste0(dir_out, "similarity_score_tissue_", today, ".rds"))



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
saveRDS(sim_score, file = paste0(dir_out, "similarity_score_tissue_cluster_EAE_", today, ".rds"))



####### UT ########
#### load data ####
treatment <- "UT"
tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL")
dir_cluster <- paste0("2_pipeline/clustering/UT_GFPall_intra/", clustering_date, "/")

## meta data
meta <- list()
for (tissue in tissue_vec) {
  tmp  <- read.table(paste0(dir_cluster, "FILES/", tissue, "/meta_data.csv"), header = T, row.names = 1, sep = ",")
  meta[[tissue]] <- tmp[,-grep("integrated_snn_res", colnames(tmp))]
}
meta_all <- Reduce(rbind, meta)

## merge clonotype info to meta data
meta_all$clonotype_id <- clty_all[rownames(meta_all),"clonotype_id"]
meta_all$pooled_clonotype_id <- paste(meta_all$treatment, meta_all$batch, meta_all$mouse, meta_all$clonotype_id, sep = "_")
meta_all$sample <- paste(meta_all$batch, meta_all$mouse, sep = "_") %>% factor()
meta_all$tissue_cluster <- factor(paste(meta_all$tissue, meta_all$seurat_clusters, sep = "_"))

meta_noNA <- meta_all[!is.na(meta_all$clonotype_id),] %>% droplevels()
stopifnot(all(meta_noNA$treatment == treatment)) ## check if treatment is correct (should be EAE cells only)

## similarity between all pairs of tissue clusters
all_clonotypes <- unique(meta_noNA$pooled_clonotype_id)
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
saveRDS(sim_score, file = paste0(dir_out, "similarity_score_tissue_cluster_UT_", today, ".rds"))

# ## compute sharing
# clones_SPL_0 <- meta_noNA$pooled_clonotype_id[meta_noNA$tissue_cluster == "SPL_0"]
# clones_SPL_1 <- meta_noNA$pooled_clonotype_id[meta_noNA$tissue_cluster == "SPL_1"]
# 
# n_shared <- meta_noNA %>% 
#   rownames_to_column(var = "cell_name") %>% 
#   mutate(in_SPL = factor((pooled_clonotype_id %in% clones_SPL_0) + 2 * (pooled_clonotype_id %in% clones_SPL_1), levels = c(0, 1, 2, 3), labels = c("none", "SPL_0", "SPL_1", "SPL_0&1"))) %>% 
#   group_by(tissue_cluster, in_SPL) %>% summarise(n_in_SPL = length(pooled_clonotype_id)) %>% ungroup() %>% 
#   mutate(tissue = gsub("_.*", "", tissue_cluster) %>% factor(levels = tissue_vec)) %>% 
#   group_by(tissue, in_SPL) %>% mutate(pct = 100 * n_in_SPL / sum(n_in_SPL)) %>% ungroup()
# 
# p_cts <- list()
# p_pct <- list()
# for (i in tissue_vec) {
#   cat("Plotting", i, "...\n")
#   tmp <- n_shared %>% filter(tissue == i) %>% droplevels()
#   fisher_out <- tmp %>% select(tissue_cluster, in_SPL, n_in_SPL) %>% 
#     spread(key = in_SPL, value = n_in_SPL, fill = 0) %>% 
#     column_to_rownames(var = "tissue_cluster") %>%
#     data.frame() %>% as.matrix() %>% 
#     fisher.test(., simulate.p.value = TRUE, B = 1e6)
#   
#   p_cts[[i]] <- ggplot(tmp) +
#     geom_bar(aes(x = in_SPL, y = n_in_SPL, fill = tissue_cluster), stat = "identity", position = "stack", width = 0.8) +
#     labs(x = "", y = "Number of cells", fill = "Cluster") +
#     ggtitle(paste0("p-value=", formatC(fisher_out$p.value, digits = 1, format = "E"))) +
#     theme_bw()
#   p_pct[[i]] <- ggplot(tmp) +
#     geom_bar(aes(x = in_SPL, y = pct, fill = tissue_cluster), stat = "identity", position = "stack", width = 0.8) +
#     labs(x = "", y = "Percentage", fill = "Cluster") +
#     theme_bw()
# }
# 
# pdf(paste0(dir_out, "share_with_SPL0_SPL1_", today, ".pdf"), width = 10, height = 5)
# for (i in tissue_vec) {
#   grid.arrange(p_cts[[i]], p_pct[[i]], nrow = 1)
# }
# dev.off()
