## configuration
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(Seurat)
library(dplyr)

today <- "2020-07-23"
dir_out <- paste0("2_pipeline/Pdl1_Arazi/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

## load data
mtx <- read.table("0_data/other/Arazi_SDY997_EXP1517/gene_by_cell_exp_mat.736297.txt", header = T) ## log reads per 10,000
meta <- read.table("0_data/other/Arazi_SDY997_EXP1517/SDY997_EXP15176_celseq_meta.tsv.725704", header = T, stringsAsFactors = F)
cell_cluster <- read.table("0_data/other/Arazi_SDY997_EXP1517/cluster_per_cell.736296.txt", header = T, stringsAsFactors = F)
stopifnot(all(cell_cluster$cell_name == colnames(mtx)))

cts <- t(t(exp(mtx) - 1) * meta$molecules[match(colnames(mtx), meta$cell_name)] / 10000) %>% round()
colnames(cts) <- rownames(meta)[match(colnames(mtx), meta$cell_name)]
meta_sub <- meta[match(cell_cluster$cell_name, meta$cell_name),]
meta_sub$cluster <- cell_cluster$cluster

so <- CreateSeuratObject(counts = cts, meta.data = meta_sub)
so <- NormalizeData(so)
so <- ScaleData(so)
so <- FindVariableFeatures(so)
so <- RunPCA(so)
ElbowPlot(so, ndims = 50)
so <- FindNeighbors(so, reduction = "pca", dims = 1:20)
so <- FindClusters(so)

so <- RunUMAP(so, reduction = "pca", dims = 1:20)

c12 <- which(so$seurat_clusters == 12)
ct6 <- which(so$cluster == "CT6")


pdf(paste0(dir_out, "cluster_overview_", today, ".pdf"))
DimPlot(so, reduction = "umap", group.by = "cluster", label = T)
DimPlot(so, reduction = "umap", group.by = "seurat_clusters", label = T)
FeaturePlot(so, features = "ISG15", reduction = 'umap')
VlnPlot(so, features = "ISG15", assay = "RNA")
FeaturePlot(so, features = "IFIT3", reduction = 'umap')
VlnPlot(so, features = "IFIT3", assay = "RNA")
FeaturePlot(so, features = "CD274", reduction = 'umap')
VlnPlot(so, features = "CD274", assay = "RNA")
vennDiagram(vennCounts(data.frame(S12 = union(c12, ct6) %in% c12, CT6 = union(c12, ct6) %in% ct6)))
dev.off()


so_CT <- so[,grepl("CT", so$cluster)]
pdf(paste0(dir_out, "cluster_overview_CT_", today, ".pdf"))
DimPlot(so_CT, reduction = "umap", group.by = "cluster", label = T)
DimPlot(so_CT, reduction = "umap", group.by = "seurat_clusters", label = T)
FeaturePlot(so_CT, features = "ISG15", reduction = 'umap')
VlnPlot(so_CT, features = "ISG15", assay = "RNA")
FeaturePlot(so_CT, features = "IFIT3", reduction = 'umap')
VlnPlot(so_CT, features = "IFIT3", assay = "RNA")
FeaturePlot(so_CT, features = "CD274", reduction = 'umap')
VlnPlot(so_CT, features = "CD274", assay = "RNA")
vennDiagram(vennCounts(data.frame(S12 = union(c12, ct6) %in% c12, CT6 = union(c12, ct6) %in% ct6)))
dev.off()

markers <- FindMarkers(so_CT, ident.1 = 12)
markers["CD274",]

cd274 <- data.frame(
  in_cluster = so_CT$seurat_clusters == 12,
  expr = so_CT[["RNA"]]@data["CD274",]
)

saveRDS(so, file = paste0(dir_out, "so_", today, ".rds"))
saveRDS(so_CT, file = paste0(dir_out, "so_CT_", today, ".rds"))
saveRDS(cd274, file = paste0(dir_out, "df_cd274_", today, ".rds"))
saveRDS(markers, file = paste0(dir_out, "markers_", today, ".rds"))


