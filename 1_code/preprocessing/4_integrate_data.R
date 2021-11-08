#### integrate all samples after QC #####

#### configuration ####
rm(list = ls())
options(future.globals.maxSize= Inf)
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-01-29"
  qc_date <- "2020-01-22"
} else {
  today <- cargs[1]
  qc_date <- cargs[2]
}
dir_out <- paste0("2_pipeline/preprocessing/integration/", today, "/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)


library(Seurat)

## read data
so_list <- readRDS(paste0("2_pipeline/preprocessing/QC/", qc_date, "/so_list.rds"))

### exclude non-hashing samples <100 cells or hashing samples <50 cells after QC
### rename cells to make sure they are unique
idx_exclude <- c()
for (i in 1:length(so_list)) {
  # cat(i,"...\n")
  if ((ncol(so_list[[i]]) < 50) | (ncol(so_list[[i]]) < 100 & so_list[[1]]$is_hashing[1])) {
    idx_exclude <- c(idx_exclude, i)
  } else {
    so_list[[i]] <- RenameCells(so_list[[i]], new.names = paste(so_list[[i]]$orig.ident, colnames(so_list[[i]]), sep = "_"))
  }
}
so_list[idx_exclude] <- NULL


### merge samples in each batch
batch_vec <- paste0("b", 1:9)
treatment_vec <- c("EAE", "UT")
mouse_vec <- c("m1", "m2")
so_list_merged <- list()
for (batch in batch_vec) {
  if (batch %in% batch_vec[1:5]) { # non-hashing batch
    sample_name <- batch
    cat("Merging", sample_name, "...\n")
    tmp <- so_list[grep(sample_name, names(so_list))]
    if (length(tmp) == 0) next
    so_list_merged[[sample_name]] <- merge(tmp[[1]], tmp[-1])
  } else { #hashing batch
    for (treatment in treatment_vec) {
      for (mouse in mouse_vec) {
        sample_name <- paste(treatment, batch, mouse, sep = "_")
        cat("Merging", sample_name, "...\n")
        tmp <- so_list[grep(sample_name, names(so_list))]
        if (length(tmp) == 0) next
        so_list_merged[[sample_name]] <- merge(tmp[[1]], tmp[-1])
      }
    }
  }
}
rm(so_list, tmp)
# saveRDS(so_list_merged, file = paste0(dir_out, "so_merged.rds"))
# so_list_merged = readRDS(file = paste0(dir_out, "so_merged.rds"))

#### integrate data from different batches
## SCTransform
so_list_merged <- lapply(so_list_merged, FUN = SCTransform, return.only.var.genes = FALSE)
saveRDS(so_list_merged, file = paste0(dir_out, "so_SCTransformed.rds"))
so_list_merged <- readRDS(file = paste0(dir_out, "so_SCTransformed.rds"))

## Integration
all_commom_features <- Reduce(intersect, lapply(so_list_merged, rownames)) ## genes
features <- SelectIntegrationFeatures(object.list = so_list_merged, nfeatures = Inf)
so_list_merged <- PrepSCTIntegration(object.list = so_list_merged,
                                     anchor.features = features,
                                     verbose = FALSE)
so_list_merged <- lapply(X = so_list_merged,
                         FUN = RunPCA,
                         verbose = FALSE,
                         features = features)
anchors <- FindIntegrationAnchors(object.list = so_list_merged,
                                      normalization.method = "SCT",
                                      reduction = "rpca",
                                      anchor.features = features,
                                      verbose = TRUE)
# saveRDS(anchors, file = paste0(dir_out, "anchors.rds"))
# anchors <- readRDS(file = paste0(dir_out, "anchors.rds"))
so_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                               verbose = TRUE)

saveRDS(so_integrated, file = paste0(dir_out, "so_integrated_SCTransformed.rds"))