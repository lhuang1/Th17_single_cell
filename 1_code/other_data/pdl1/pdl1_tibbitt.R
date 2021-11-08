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


today <- "2020-07-23"
dir_out <- paste0("2_pipeline/Pdl1_Tibbitt/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

d1 <- read.table("0_data/other/Tibbitt_GSE131935_RAW/GSE131935_SS2_15_0160_rpkms.tab.gz", header = T, stringsAsFactors = F)
colnames(d1) <- paste0("SS2_15_0160_star_mm10_", colnames(d1))
d2 <- read.table("0_data/other/Tibbitt_GSE131935_RAW/GSE131935_SS2_17_218_rpkms_.tab.gz", header = T, stringsAsFactors = F)
colnames(d2) <- paste0("SS2_17_218_star_mm10_", colnames(d2))
d3 <- read.table("0_data/other/Tibbitt_GSE131935_RAW/GSE131935_SS2_17_449_rpkms.tab.gz", header = T, stringsAsFactors = F)
colnames(d3) <- paste0("SS2_17_449_star_mm10_", colnames(d3))

stopifnot(all(d1[,1] == d2[,1]) & all(d1[,1] == d3[,1]))

d <- cbind(d1, cbind(d2[,-1], d3[,-1]))
colnames(d)[1] <- "gene"
# rm(d1, d2, d3)

# ## check expression for duplicated gene names
dup_idx <- duplicated(d$gene) %>% which()
dup_gene <- unique(d$gene[dup_idx])
# dup_dist <- sapply(dup_gene, function(i) {
#   d_sub <- d[d$gene == i,-1]
#   dist(d_sub) %>% as.numeric() %>% sum
# })
# dup_total <- sapply(dup_gene, function(i) {
#   d_sub <- d[d$gene == i,-1]
#   sum(d_sub)
# })
# all(which(dup_dist > 0) == which(dup_total > 0))
# # seems that if a gene is detected, then two rows can be different, thus take the summation
## add up duplicated genes
for (g in dup_gene) {
  d[d$gene == g, -1] <- rep(colSums(d[d$gene == g, -1]), sum(d$gene == g)) %>% matrix(., ncol = ncol(d) - 1, byrow = T)
}
d <- d[-dup_idx,]
rownames(d) <- d$gene
d <- d[,colnames(d) != "gene"]
d <- as.sparse(d)
saveRDS(d, file = paste0(dir_out, "count_matrix_", today, ".rds"))

# d <- readRDS(file = paste0(dir_out, "count_matrix_", today, ".rds"))


#### preprocessing
## load sample info
sample_info <- read.table("0_data/other/Tibbitt_GSE131935_RAW/GSE131935_series_matrix.txt.gz", sep = "\t", skip = 30, nrows = 8, header = T, row.names = 1, stringsAsFactors = F)
sample_info <- t(sample_info) %>% data.frame() %>% 
  filter(X.Sample_source_name_ch1 %in% c("Day 15 BAL T helper cells mouse 1", "Day 15 BAL T helper cells mouse 2"))
head(sample_info)


#### QC
## check batch effects
d <- d[,rownames(sample_info)]
so <- CreateSeuratObject(counts = d, 
                         meta.data = sample_info)
so$mouse <- sub("Day 15 BAL T helper cells ", "", so$X.Sample_source_name_ch1)
## normalization & clustering
so <- NormalizeData(so, assay = "RNA")
so <- ScaleData(so)
so <- FindVariableFeatures(so)
so <- RunPCA(so)
ElbowPlot(so, ndims = 50)
DimPlot(so, reduction = "pca", group.by = "mouse")
so <- FindNeighbors(so, reduction = "pca", dims = 1:20)
so <- FindClusters(so)

so <- RunUMAP(so, reduction = "pca", dims = 1:20)

pdf(paste0(dir_out, "cluster_overview_", today, ".pdf"))
DimPlot(so, reduction = "umap", group.by = "mouse")
DimPlot(so, reduction = "umap", group.by = "seurat_clusters", label = T)
FeaturePlot(so, features = "Isg15", reduction = 'umap')
VlnPlot(so, features = "Isg15", assay = "RNA")
FeaturePlot(so, features = "Ifit3", reduction = 'umap')
VlnPlot(so, features = "Ifit3", assay = "RNA")
FeaturePlot(so, features = "Cd274", reduction = 'umap')
VlnPlot(so, features = "Cd274", assay = "RNA")
dev.off()


markers <- FindMarkers(so, ident.1 = 4)

cd274 <- data.frame(
  in_cluster = so$seurat_clusters == 4,
  expr = so[["RNA"]]@data["Cd274",]
)

saveRDS(so, file = paste0(dir_out, "so_", today, ".rds"))
saveRDS(cd274, file = paste0(dir_out, "df_cd274_", today, ".rds"))
saveRDS(markers, file = paste0(dir_out, "markers_", today, ".rds"))

