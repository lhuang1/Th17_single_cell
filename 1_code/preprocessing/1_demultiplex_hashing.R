#### parse demuxEM output for cell-hashing batches ####
#### configuration ####
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
today <- Sys.Date()
convert_date <- "2020-01-22"
dir_out <- paste0("2_pipeline/preprocessing/demultiplex/", today, "/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)


library(Seurat)
library(ggplot2)
source("1_code/utils.R")
source("1_code/preprocessing/utils.R")


library(reticulate)
ad <- import("anndata", convert = FALSE)
pd <- import("pandas", convert = FALSE)


batch_vec <- paste0("batch", 6:9)
mouse_vec <- c("m1", "m2")
treatment_vec <- c("UT", "EAE")

for (batch in batch_vec) {
  for (treatment in treatment_vec) {
    for (mouse in mouse_vec) {
      sample_name <- paste(batch, "hashed", treatment, mouse, sep = "_")
      dir_hto <- paste0(dir_proj, "data/single_cell/demuxEM/", sample_name, "/", treatment, "_", mouse, "_ADTs.h5ad")
      dir_rna <- paste0(dir_proj, "data/single_cell/demuxEM/", sample_name, "/", treatment, "_", mouse, "_demux.h5ad")
      if ((!file.exists(dir_hto)) | !file.exists(dir_rna)) next
      cat("Processing", sample_name, "...\n")
      # hashtag UMI counts
      hto <- ReadH5AD(dir_hto, assay = "HTO")
      # hashtag assignment
      rna_adata <- ad$read_h5ad(dir_rna)
      rna_meta <- py_to_r(rna_adata$obs)
      # only keep valid cells (hto contains hundreds of thousands of "cells")
      common_cells <- intersect(rownames(rna_meta), colnames(hto))
      hto <- hto[,common_cells] 
      hto$assignment_hashtag <- as.character(rna_meta[common_cells, "assignment"])
      rm(rna_adata, rna_meta)
      
      ## parse DemuxEM results
      hto$assignment_tissue <- parse_assignment(hto$assignment_hashtag, batch = batch)
      
      ## dimension reduction and visualization
      hto <- NormalizeData(hto, normalization.method = "CLR")
      hto <- RunUMAP(hto, features = rownames(hto))
      
      p1 <- ggplot() +
        geom_bar(data = hto@meta.data, aes(x = assignment_tissue, fill = assignment_tissue)) + 
        geom_text(data = as.data.frame(table(hto$assignment_tissue)),
                  aes(x = Var1, y = Freq/2, label = Freq)) +
        labs(x = "", y = "Number of cells", fill = "Assignment") +
        theme_bw()
      
      p2 <- UMAPPlot(hto, group.by = "assignment_tissue") + theme_bw()
      
      p3 <- FeaturePlot(hto, features = rownames(hto), reduction = 'umap', pt.size = 0.2)
      
      pdf(paste0(dir_out, "UMAP_demuxEM_", sample_name, ".pdf"))
      print(p1)
      print(p2)
      print(p3)
      dev.off()
      
      saveRDS(hto@meta.data, file = paste0(dir_out, "meta_", sample_name, ".rds"))
      write.csv(hto@meta.data, file = paste0(dir_out, "meta_", sample_name, ".csv"), row.names = T, quote = F)
    }
  }
}
