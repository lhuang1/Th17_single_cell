#### TCR sharing between SPL0 and SPL1
rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(dplyr)
library(tibble)
library(tidyr)

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()
  clustering_date <- "2020-03-25"
  clonotyping_type_date <- "dominant_2020-04-02"
} else {
  today <- cargs[1]
  clustering_date <- cargs[2]
  clonotyping_date <- cargs[3]
}

dir_out <- paste0("2_pipeline/TCR/EAE_SPL0_vs_SPL1/")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

#### load data ####
dir_spl_cluster <- paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/")
meta <- read.table(paste0(dir_spl_cluster, "FILES/meta_data_switched.txt"), header = T, sep = " ")

## clonotypes
clty_all <- read.csv(paste0("2_pipeline/TCR/clonotype_assignment_", clonotyping_type_date, ".csv"), row.names = 1)

## merge clonotype info to meta data
meta$clonotype_id <- clty_all[rownames(meta),"clonotype_id"]
meta$pooled_clonotype_id <- paste(meta$treatment, meta$batch, meta$mouse, meta$clonotype_id, sep = "_")
meta$tissue_cluster <- factor(paste(meta$tissue, meta$seurat_clusters, sep = "_"))

meta_noNA <- meta[!is.na(meta$clonotype_id),] %>% droplevels()

clones_SPL_0 <- meta_noNA$pooled_clonotype_id[meta_noNA$tissue_cluster == "SPL_0"] %>% unique()
clones_SPL_1 <- meta_noNA$pooled_clonotype_id[meta_noNA$tissue_cluster == "SPL_1"] %>% unique()
clones_SPL_0_and_1 <- intersect(clones_SPL_0, clones_SPL_1)
clones_SPL_0_or_1 <- union(clones_SPL_0, clones_SPL_1)

col_spl_0_1 <- SPL_EAE_cluster_0_1_colors()


## clone size in SPL_0 and SPL_1, stratify by shared or not
size_df <- meta_noNA %>% filter(tissue_cluster %in% c("SPL_0", "SPL_1")) %>%
  group_by(pooled_clonotype_id, tissue_cluster) %>% tally(name = "clone_size") %>%
  arrange(clone_size) %>%
  mutate(type = ifelse(pooled_clonotype_id %in% clones_SPL_0_and_1, "Shared", "Unique")) %>%
  mutate(type = factor(type, levels = c("Shared", "Unique")),
         clone_size_cat = factor(cut(clone_size, c(0:5, Inf), labels = c(1:5, "6+")), levels = rev(c(1:5, "6+"))))
## pct of SPL_0 an SPL_1 cells sharing TCR 
pct_df <- size_df %>% 
  group_by(tissue_cluster) %>% 
  summarise(total_cells = sum(clone_size), 
            pct_shared = sum(clone_size[type == "Shared"]) / total_cells * 100,
            pct_not_shared = 100 - pct_shared) %>% 
  gather(key = "is_shared", value = "pct", -tissue_cluster, -total_cells) %>% 
  mutate(is_shared = factor(is_shared, levels = c("pct_not_shared", "pct_shared"), labels = c("Unique", "Shared")))

saveRDS(size_df, file = paste0(dir_out, "clone_size_SPL0_SPL1_shared_unique_", today, ".rds"))
saveRDS(pct_df, file = paste0(dir_out, "cell_pct_SPL0_SPL1_shared_unique_", today, ".rds"))
