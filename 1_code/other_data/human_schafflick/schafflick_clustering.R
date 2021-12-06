####### Analyze human (Schafflick) data ########

## configuration
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(dplyr)
library(tibble)
library(openxlsx)
source("1_code/utils.R")
source("1_code/clustering/clustering_utils.R")

integrate_date <- "2020-06-01"
today <- "2020-06-11"
dir_out <- paste0("2_pipeline/human_shafflick/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

##### load data ######
pbmc <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/so_combined_PBMCs.rds"))
pbmc@meta.data <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/meta_combined_clustered_PBMCs.rds"))
pbmc@reductions <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/reductions_combined_clustered_PBMCs.rds"))
csf <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/so_combined_CSF.rds"))
csf@meta.data <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/meta_combined_clustered_CSF.rds"))
csf@reductions <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/reductions_combined_clustered_CSF.rds"))


#### check how treatment/individuals mix after integration
pdf(paste0("2_pipeline/human_shafflick/", integrate_date, "/UMAP_mixing.pdf"))
DimPlot(csf, group.by = "is_MS", label = F) + ggtitle("CSF")
DimPlot(csf, group.by = "individual", label = F) + ggtitle("CSF")
DimPlot(pbmc, group.by = "is_MS", label = F) + ggtitle("PBMC")
DimPlot(pbmc, group.by = "individual", label = F) + ggtitle("PBMC")
dev.off()


# restrict to CD4+ cells (excluding Tregs)
pbmc <- pbmc[,pbmc$seurat_clusters %in% c(0, 1, 4, 10)] ## c(1, 0, 7)


csf <- csf[,csf$seurat_clusters %in% c(0, 1, 2)] ## c(0, 3)



##### clustering and dimension reduction on CD4 ######
pbmc <- RunPCA(pbmc, npcs = 15, verbose = FALSE)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:15, min.dist = 0.5)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)
saveRDS(pbmc@meta.data, file = paste0(dir_out, "meta_CD4_PBMCs.rds"))
saveRDS(pbmc@reductions, file = paste0(dir_out, "reductions_CD4_PBMCs.rds"))
pbmc_markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.3, max.cells.per.ident = 500,
                               logfc.threshold = 0.25)
saveRDS(pbmc_markers, paste0(dir_out, "pbmc_markers.rds"))
save_and_plot_markers(pbmc, pbmc_markers, dir_markers = dir_out, tissue = "PBMCs", pt.size = 0.5)


csf <- RunPCA(csf, npcs = 15, verbose = FALSE)
csf <- RunUMAP(csf, reduction = "pca", dims = 1:15, min.dist = 0.5)
csf <- FindNeighbors(csf, reduction = "pca", dims = 1:15)
csf <- FindClusters(csf, resolution = 0.5)
saveRDS(csf@meta.data, file = paste0(dir_out, "meta_CD4_CSF.rds"))
saveRDS(csf@reductions, file = paste0(dir_out, "reductions_CD4_CSF.rds"))
csf_markers <- FindAllMarkers(object = csf, only.pos = TRUE, min.pct = 0.3, max.cells.per.ident = 500,
                               logfc.threshold = 0.25)
saveRDS(csf_markers, paste0(dir_out, "csf_markers.rds"))
save_and_plot_markers(csf, csf_markers, dir_markers = dir_out, tissue = "CSF", pt.size = 0.5)

pdf(paste0("2_pipeline/human_shafflick/", today, "/UMAP_clusters_CD4.pdf"))
DimPlot(csf, group.by = "seurat_clusters", label = T) + ggtitle("CSF")
DimPlot(pbmc, group.by = "seurat_clusters", label = T) + ggtitle("PBMC")
dev.off()


#### check again in CD4 cells how treatment/individuals mix after integration
pdf(paste0("2_pipeline/human_shafflick/", today, "/UMAP_mixing.pdf"))
DimPlot(csf, group.by = "is_MS", label = F) + ggtitle("CSF")
DimPlot(csf, group.by = "individual", label = F) + ggtitle("CSF")
DimPlot(pbmc, group.by = "is_MS", label = F) + ggtitle("PBMC")
DimPlot(pbmc, group.by = "individual", label = F) + ggtitle("PBMC")
dev.off()

