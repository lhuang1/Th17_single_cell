####### Analyze human (Schafflick) data ########

## configuration
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggridges)
library(dplyr)
library(tibble)
library(openxlsx)
source("1_code/utils.R")

integrate_date <- "2020-06-01"
clustering_date <- "2020-06-11"
today <- "2020-07-01"
dir_out <- paste0("2_pipeline/human_shafflick/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

##### load data ######
pbmc <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/so_combined_PBMCs.rds"))
meta <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/meta_CD4_PBMCs.rds"))
pbmc <- pbmc[,rownames(meta)]
pbmc@meta.data <- meta
pbmc@reductions <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/reductions_CD4_PBMCs.rds"))
DefaultAssay(pbmc) <- "integrated"
csf <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/so_combined_CSF.rds"))
meta <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/meta_CD4_CSF.rds"))
csf <- csf[,rownames(meta)]
csf@meta.data <- meta
csf@reductions <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/reductions_CD4_CSF.rds"))
DefaultAssay(csf) <- "integrated"
## take out non-CD4  
if (integrate_date == 2020-06-01 & clustering_date == 2020-06-11) {
  pbmc <- pbmc[,!(pbmc$seurat_clusters %in% c(2, 5, 7, 8))] 
  csf <- csf[,!(csf$seurat_clusters %in% c(2, 4))]
}

##### Signatures ######
## load signature data
# Gaublomme's
sigs_raw_gaublomme <- read.table("0_data/gene_lists/all_Gaublomme_2016_signatures.pgf.txt", stringsAsFactors = F)
sigs_raw_gaublomme$sig_name <- paste0(sigs_raw_gaublomme$V1, "-", sigs_raw_gaublomme$V2)
sigs_gaublomme <- lapply(unique(sigs_raw_gaublomme$sig_name), function(i){sigs_raw_gaublomme$V3[sigs_raw_gaublomme$sig_name == i]})
names(sigs_gaublomme) <- unique(sigs_raw_gaublomme$sig_name)
# CXCR6 vs. SLAMF6
bulk_sigs <- readRDS("2_pipeline/spl_clusters/bulk/2020-02-04/signature_slamf6_vs_cxcr6_all.rds")
bulk_sigs <- mouse_to_human(bulk_sigs)
sc_sigs <- readRDS("2_pipeline/differential_expression/SPL_EAE_0_1/2020-04-01/sc_signature.rds")
sc_sigs <- mouse_to_human(sc_sigs)
# proliferating from our clustering
sigs_raw_proliferating <- readRDS("2_pipeline/differential_expression/Cluster_UT_GFPpos_inter/2020-04-16/results.rds")
sigs_proliferating <- sigs_raw_proliferating[["7"]] %>%
  rownames_to_column(var = "gene") %>% filter(logFC > log2(1.5) & FDR < 0.05) %>% select(gene) %>% unlist %>% as.character() %>%
  mouse_to_human()

sigs_other <- list()
sigs_other[['Cxcr6_bulk-plus']] <- bulk_sigs$"Cxcr6-plus"
sigs_other[['Slamf6_bulk-plus']] <- bulk_sigs$"Slamf6-plus"
sigs_other[['Cxcr6_sc-plus']] <- sc_sigs$"Cxcr6_sc-plus"
sigs_other[['Slamf6_sc-plus']] <- sc_sigs$"Slamf6_sc-plus"
sigs_other[['Pathogenic_Th17']] <- c("Cxcl3", "Il22", "Il3", "Ccl4", "Gzmb", "Lrmp", "Ccl5", "Casp1", "Csf2", "Ccl3", "Tbx21", "lcos", "ll7r", "Stat4", "Lgals3", "Lag3") %>% mouse_to_human()
sigs_other[["Proliferating"]] <- sigs_proliferating
signature_list <- c(sigs_other, sigs_gaublomme)
saveRDS(signature_list, file = paste0(dir_out, "signature_list.rds"))

#### compute signature expression score
sapply(signature_list, function(x){
  length(intersect(x, rownames(pbmc)))
})
pbmc <- AddModuleScore(pbmc, features = signature_list, search = F)
colnames(pbmc@meta.data)[grep("Cluster", colnames(pbmc@meta.data))] <- names(signature_list)
pbmc$tissue_MS <- paste(pbmc$tissue, pbmc$is_MS, sep = ".")

sapply(signature_list, function(x){
  length(intersect(x, rownames(csf)))
})
csf <- AddModuleScore(csf, features = signature_list, search = F)
colnames(csf@meta.data)[grep("Cluster", colnames(csf@meta.data))] <- names(signature_list)
csf$tissue_MS <- paste(csf$tissue, csf$is_MS, sep = ".")

saveRDS(pbmc@meta.data, file = paste0(dir_out, "meta_pbmc.rds"))
saveRDS(csf@meta.data, file = paste0(dir_out, "meta_csf.rds"))
    
