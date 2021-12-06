#################################################
#### Assign clonotype to cells with TCR info ####

## perfect match of CDR3 nucleotide sequences
## assign for each mouse independently

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(dplyr)
library(tibble)

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()
  combine_tcr_date <- "2020-03-24"
} else {
  today <- cargs[1]
  combine_tcr_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/TCR/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

### load data
meta_dominant <- readRDS(paste0("2_pipeline/preprocessing/meta_processed_dominant_TCR_", combine_tcr_date, ".rds"))

### assign clonotype
treatment_vec <- c("UT", "EAE")
batch_vec <- c("b6", "b7", "b8", "b9")
mouse_vec <- c("m1", "m2")
meta_clonotype <- data.frame(matrix(NA, nrow = sum(!is.na(meta$raw_clonotype_id)), ncol = 2))
rownames(meta_clonotype) <- rownames(meta)[!is.na(meta$raw_clonotype_id)]
colnames(meta_clonotype) <- c("clonotype_id", "clone_size")
for (treatment in treatment_vec) {
  for (batch in batch_vec) {
    for (mouse in mouse_vec) {
      meta_sub <- meta[meta$treatment == treatment & meta$batch == batch & meta$mouse == mouse & !is.na(meta$raw_clonotype_id),]
      if (nrow(meta_sub) == 0) next ## no such treatment, batch, mouse combination
      cat(treatment, batch, mouse, "\n")
      clonotypes <- meta_sub %>% rownames_to_column(var = "cellname") %>% 
        group_by(TRA, TRB) %>% summarize(cells_in_clone = paste(unique(cellname), collapse = ";"), clone_size = length(unique(cellname))) %>% 
        arrange(-clone_size) %>%
        data.frame()
      clonotypes$clonotype_id <- paste("clonotype", seq_along(clonotypes[,1]), sep = "")
      cells_list <- strsplit(clonotypes$cells_in_clone, split = ";")
      for (i in seq_along(clonotypes$clonotype_id)) {
        meta_clonotype[cells_list[[i]],"clonotype_id"] <-  clonotypes$clonotype_id[i]
        meta_clonotype[cells_list[[i]],"clone_size"] <-  clonotypes$clone_size[i]
      }
    }
  }
}

write.csv(meta_clonotype, file = paste0(dir_out, "clonotype_assignment_dominant_", today, ".csv"))

