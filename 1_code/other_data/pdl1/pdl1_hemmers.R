rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")
library(Seurat)
library(dplyr)


today <- "2020-07-23"
dir_out <- paste0("2_pipeline/Pdl1_Hemmers/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

## load data
mtx <- read.csv("0_data/other/Hemmers_GSE134902_RAW/GSM3978656_MatureCD4SP_dense.csv.gz",
                header = T, row.names = 1)

expr <- mtx[,-1] %>% t
meta <- data.frame(
  Cluster = mtx[,1],
  row.names = rownames(mtx)
)
so <- CreateSeuratObject(counts = expr,
                         meta.data = meta)
# rm(mtx)

so <- so[,so$nCount_RNA >= 1000]

so <- SCTransform(so)
so <- RunPCA(so, assay = "SCT")
ElbowPlot(so, ndims = 50)
so <- FindNeighbors(so, reduction = "pca", dims = 1:20)
so <- FindClusters(so)

so <- RunUMAP(so, reduction = "pca", dims = 1:15)
DimPlot(so, label = T)
FeaturePlot(so, features = "ISG15")
FeaturePlot(so, features = "CD274")

markers <- FindMarkers(so, ident.1 = 5)


so <- NormalizeData(so, assay = "RNA")

cd274 <- data.frame(
  in_cluster = so$seurat_clusters == 5,
  expr = so[["RNA"]]@data["CD274",]
)

saveRDS(so, file = paste0(dir_out, "so_", today, ".rds"))
saveRDS(cd274, file = paste0(dir_out, "df_cd274_", today, ".rds"))
saveRDS(markers, file = paste0(dir_out, "markers_", today, ".rds"))