library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)

dir_proj <- "/singerlab/linglin/Th17_single_cell_eae_ut/"
clustering_date <- "2020-03-25"
prep_type_date <- "dominant_TCR_2020-03-24"
meta <- readRDS(paste0(dir_proj, "2_pipeline/preprocessing/meta_processed_", prep_type_date, ".rds"))
## load clonotype data
clonotyping_type_date <- "dominant_2020-04-02"
clty_all <- read.csv(paste0(dir_proj, "2_pipeline/TCR/clonotype_assignment_", clonotyping_type_date, ".csv"), 
                     row.names = 1)
## merge clonotypes IDs to meta data
meta$clonotype_id <- clty_all[rownames(meta),"clonotype_id"]
meta_noNA <- meta[!is.na(meta$clonotype_id),] %>% droplevels()
meta_noNA$tissue_cluster <- paste(meta_noNA$tissue, meta_noNA$seurat_clusters, sep = "_")
meta_noNA$pooled_cid <- paste(meta_noNA$batch, meta_noNA$mouse, meta_noNA$clonotype_id, sep = "_")

## compare mice
## UT
meta_noNA %>% filter(treatment == "UT") %>% 
  group_by(batch, mouse) %>% summarise(n_clones = length(unique(clonotype_id)))
# 774 + 897 + 1014 + 1033 = 3718
# 3718 - 11 - 1 = 3706

## EAE
meta_noNA %>% filter(treatment == "EAE") %>% 
  group_by(batch, mouse) %>% summarise(n_clones = length(unique(clonotype_id)))
# 2085 + 1639 + 415 +1151 = 5290
# 5290 - 12 = 5278

# 3718 + 5290 == 9008
# 9008 - 37 - 1 = 8970 (-1 because 1 clonotype is shared by 3 mice)

meta_noNA %>% filter(treatment == "EAE") %>% 
  group_by(batch, mouse) %>% summarise(frac_share =mean(rowSums(table(clonotype_id, tissue)) > 1))
meta_sub <- meta_noNA %>% filter(treatment == "EAE") 
mean(rowSums(table(meta_sub$pooled_cid, meta_sub$tissue)) > 1)

# based on table
# UT-UT: 11
# EAE-EAE: 12
# UT-EAE: 13
# UT-UT-EAE: 1
# total: 37

