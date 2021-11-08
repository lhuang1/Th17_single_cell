#### compute Il23 related signature scores
rm(list = setdiff(ls(), "so_all"))

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(gridExtra)
library(ggsignif)
library(openxlsx)
set.seed(1)
source("1_code/utils.R")

#### configuration ####
today <- "2020-11-06"
clustering_date <- "2020-03-25"
prep_type_date <- "dominant_TCR_2020-03-24"


dir_out <- paste0("2_pipeline/spl_clusters/Il23/")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

#### load data ####
# so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_", prep_type_date, ".rds"))
meta <- read.table(paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/FILES/meta_data.txt"), header = T, sep = " ")

so <- so_all[,rownames(meta)[meta$seurat_clusters %in% c(0, 1)]]
so@meta.data <- meta[meta$seurat_clusters %in% c(0, 1),] %>% droplevels()
so$tissue_cluster <- paste0("SPL_", so$seurat_clusters)


## read and parse signatures
sp <- read.table("2_pipeline/COMPASS/other/signatures_from_AW/supervised_partition.csv",
                 sep = ",", stringsAsFactors = F, header = T)
sp$symbol_new <- rownames(so)[match(toupper(sp$symbol), toupper(rownames(so)))]
sp <- sp[!is.na(sp$symbol_new),]
sigs_to_plot_allon <- c("23R_KO_B623_1", "23R_KO_IL-12_1", "IL23", "IL6_IL1_IL23_1", "23R_KO_IL-12_IL-23_1", "SGK1-IL23")
sigs_parsed_allon <- lapply(sigs_to_plot_allon, function(x){
  sp$symbol_new[sp$sig_name == x & sp$direction == 1]
})
names(sigs_parsed_allon) <- sigs_to_plot_allon

fname <- "0_data/gene_lists/signatures/IL23_signatures.xlsx"
sigs_to_plot_alex <- getSheetNames(fname)
sigs_parsed_alex <- lapply(sigs_to_plot_alex, function(x){
  tmp <- read.xlsx(fname, sheet = x, rowNames = F, colNames = F)[,1]
  if (grepl("human", x)) {
    genes <- human_to_mouse(tmp)
  } else {
    genes <- rownames(so)[match(toupper(tmp), toupper(rownames(so)))] %>% setdiff(., NA)
  }
  genes
})
names(sigs_parsed_alex) <- sigs_to_plot_alex

sigs_parsed <- c(sigs_parsed_alex, sigs_parsed_allon)
saveRDS(sigs_parsed, file = "0_data/gene_lists/signatures/Il23_signatures_parsed.rds")

## compute signature scores and make plots
so <- AddModuleScore(so, sigs_parsed, name = "SIG")
sig_score_df <- so@meta.data[, grep("SIG", colnames(so@meta.data))]
colnames(sig_score_df) <- names(sigs_parsed)

saveRDS(sig_score_df, file = paste0(dir_out, "Il23_signature_scores_", today, ".rds"))


