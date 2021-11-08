#### cell and gene quality control #####
## Input: cellranger count filtered feature matrix
## Output: one big seurat object of all samples
## 1. cell QC: filter out cells that
##  (1) have non-zero expression in <500 genes
##  (2) have >5% mitochondrial genes
##  (3) have non-zero expression in <60% of the house keeping genes
## 2. gene QC: filter out genes that are not detected in any cell
## 3. merge GFP: 
##  (1) add in GFP counts
##  (2) mark cell status according to GFP +/-


#### configuration ####
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
today <- "2020-01-22" # Sys.Date()
demux_date <- "2020-01-22"
dir_out <- paste0("2_pipeline/preprocessing/QC/", today, "/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)


library(Seurat)
library(ggplot2)
source("1_code/utils.R")
source("1_code/preprocessing/utils.R")

#### QC ####
## read in house keeping genes
hk_genes <- human_to_mouse(read.table("0_data/gene_lists/HK_Satija.txt", stringsAsFactors = F)[,1])
tissue_vec <- c("SPL", "PP", "MLN", "SI", "COL", "CNS", "DLN")
treatment_vec <- c("UT", "EAE")
batch_vec <- paste0("batch", 1:9)
mouse_vec <- c("m1", "m2")

so_list <- list()
for (batch in batch_vec) {
  for (treatment in treatment_vec) {
    if (batch %in% paste0("b", 1:5)) { ## non-hashing batches
      for (tissue in tissue_vec) {
        project_name <- paste(tissue, treatment, batch, sep = "_")
        cat(project_name, " started...\n")
        so_list[[project_name]] <- prep_data_nh(tissue, treatment, batch, project_name, dir_out)
        cat(project_name, " ended!\n")
      }
    } else if (batch %in% paste0("b", 6:8)) {
      for (mouse in mouse_vec) {
        project_name <- paste(treatment, batch, mouse, sep = "_")
        cat(project_name, " started...\n")
        so_by_tissue <- prep_data_h(treatment, batch, mouse, project_name, dir_out, demux_date)
        so_list <- c(so_list, so_by_tissue)
        cat(project_name, " ended!\n")
      }
    }
    
  }
}
so_list[sapply(so_list, is.null)] <- NULL
saveRDS(so_list, file = paste0(dir_out, "so_list.rds"))

