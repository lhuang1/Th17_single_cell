#################################################
#### Assign clonotype to cells with TCR info ####

## perfect match of CDR3 nucleotide sequences
## assign for each mouse independently

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tibble)
set.seed(1)
source("1_code/utils.R")
source("1_code/tcr/utils.R")

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-04-02"
  combine_tcr_date <- "2020-03-24"
} else {
  today <- cargs[1]
  combine_tcr_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/TCR/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

### load data
meta_dominant <- readRDS(paste0("2_pipeline/preprocessing/meta_processed_dominant_TCR_", combine_tcr_date, ".rds"))
meta_2 <- readRDS(paste0("2_pipeline/preprocessing/meta_processed_withTCR_", combine_tcr_date, ".rds"))

meta_cns = read.csv("2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/CNS/meta_data.csv", header = T, row.names = 1)
ggplot(meta_cns) + geom_boxplot(aes(x = factor(seurat_clusters),  y = nFeature_RNA)) 
### assign clonotype
clonotypes_d <- meta_to_clonotype(meta_dominant)
clonotypes_2 <- meta_to_clonotype(meta_2)
clonotypes_r <- meta_dominant %>% rownames_to_column(var = "cellname") %>% 
  filter(!is.na(raw_clonotype_id)) %>% 
  group_by(orig.ident, raw_clonotype_id) %>% add_tally(name = "clone_size") %>% 
  select(raw_clonotype_id, clone_size, cellname) %>% 
  column_to_rownames(var = "cellname") %>% 
  mutate(type = "raw", clonotype_id = raw_clonotype_id) %>% 
  select(clonotype_id, clone_size, type)
  
plt_df <- rbind(rbind(cbind(clonotypes_d, type = "dominant"), cbind(clonotypes_2, type = "2ab")), clonotypes_r)
ggplot(plt_df) +
  geom_histogram(aes(x = clone_size)) +
  facet_wrap(~type) +
  theme_bw()

write.csv(clonotypes_d, file = paste0(dir_out, "clonotype_assignment_dominant_", today, ".csv"))
write.csv(clonotypes_2, file = paste0(dir_out, "clonotype_assignment_2ab_", today, ".csv"))

