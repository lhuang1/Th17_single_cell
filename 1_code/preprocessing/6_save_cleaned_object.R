### Generate a Seurat object of clean cells for downstream analysis
# also add TCR info to meta.data.

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
set.seed(1)

#### configuration ####
today <- "2020-03-24"
integration_date <- "2020-01-29"
expr_cleaning_date <- "2020-01-29"
tcr_cleaning_date <- "2020-03-24"

dir_out <- paste0("2_pipeline/preprocessing/")
if (!dir.exists(dir_out)) stop("dir_out does not exist")

### load data
so_all <- readRDS(paste0("2_pipeline/preprocessing/integration/", integration_date, "/so_integrated_SCTransformed.rds"))
keep_cells_clean <- readRDS(paste0("2_pipeline/preprocessing/cleaning/", expr_cleaning_date, "/round_3/keep_cell_names.rds"))
so_all <- so_all[,keep_cells_clean]
so_all$mouse <- sapply(strsplit(as.character(so_all$orig.ident), split = "_"), function(x){x[4]}) %>% factor()


batch_vec <- paste0("b", 6:9)
treatment_vec <- c("EAE", "UT")
mouse_vec <- c("m1", "m2")

all_contigs <- list()
for (batch in batch_vec) {
  for (treatment in treatment_vec) {
    for (mouse in mouse_vec) {
      sample_name <- paste0(treatment, "_", batch, "_", mouse)
      f_vdj <- paste0("2_pipeline/preprocessing/tcr/", tcr_cleaning_date, "/", sample_name, "_filtered_contigs_dominant.rds")
      if (!file.exists(f_vdj)) next()
      all_contigs[[sample_name]] <- readRDS(f_vdj) %>% mutate(sample_name = sample_name)
    }
  }
}
all_contigs <- Reduce(rbind, lapply(all_contigs, data.frame))

cells_expr <- sub(".*?_", "", colnames(so_all)[so_all$is_hashing])
cells_tcr <- paste(all_contigs$sample_name, all_contigs$barcode_nosuffix, sep = "_")

mean(cells_expr %in% cells_tcr) # percent have TCR (after dominant chain filtering)


tcr_df <- sapply(seq_along(cells_expr), function(i) {
  x <- cells_expr[i]
  if (x %in% cells_tcr) {
    tmp <- data.frame(all_contigs[cells_tcr == x,])
    stopifnot(length(unique(tmp$raw_clonotype_id)) == 1)
    TRA <- paste(as.character(tmp[tmp$chain == "TRA", "cdr3_nt"]), collapse = ";")
    TRB <- paste(as.character(tmp[tmp$chain == "TRB", "cdr3_nt"]), collapse = ";")
    nTRA <- sum(tmp$chain == "TRA")
    nTRB <- sum(tmp$chain == "TRB")
    raw_clonotype_id <- as.character(tmp$raw_clonotype_id[1])
    return(c(TRA, TRB, nTRA, nTRB, raw_clonotype_id))
  } else {
    return(c(NA, NA, NA, NA, NA))
  }
}) %>% t()
rownames(tcr_df) <- colnames(so_all)[so_all$is_hashing]

### merge to meta.data
so_all$TRA <- NA
so_all$TRB <- NA
so_all$nTRA <- NA
so_all$nTRB <- NA
so_all$raw_clonotype_id <- NA
so_all@meta.data[rownames(tcr_df), c("TRA", "TRB", "nTRA", "nTRB", "raw_clonotype_id")] <- tcr_df
so_all$nTRA <- as.numeric(so_all$nTRA)
so_all$nTRB <- as.numeric(so_all$nTRB)

timestamp()
cat("save so_all: ")
saveRDS(so_all, file = paste0(dir_out, "so_processed_dominant_TCR_", today, ".rds"))
saveRDS(so_all@meta.data, file = paste0(dir_out, "meta_processed_dominant_TCR_", today, ".rds"))
timestamp()
