##### process data, clustering and stain for Cxcl16 #####

### configuration
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(Seurat)
library(ggplot2)
library(dplyr)
library(openxlsx)

today <- "2020-02-26"
dir_out <- paste0("2_pipeline/cxcl16/", today, "/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

### load data
sample_info <- c(
"GSM3717013", "Priming_replicate1",
"GSM3717014", "Priming_replicate2",
"GSM3717015", "Priming_replicate3",
"GSM3717016", "CFA_replicate1",
"GSM3717017", "Acute_replicate1",
"GSM3717018", "Acute_replicate2",
"GSM3717019", "Acute_replicate3",
"GSM3717020", "Naive_replicate1",
"GSM3717021", "Naive_replicate2",
"GSM3717022", "Remitting_replicate1",
"GSM3717023", "Remitting_replicate2",
"GSM3717024", "Remitting_replicate3",
"GSM3717025", "CFA_replicate2",
"GSM3717026", "CFA_replicate3",
"GSM3717027", "Naive_replicate3",
"GSM4002707", "Acute_replicate4",
"GSM4002708", "Acute_replicate5",
"GSM4002709", "Acute_replicate6",
"GSM4002710", "Naive_replicate4",
"GSM4002711", "Remitting_replicate4",
"GSM4002712", "Remitting_replicate5",
"GSM4002713", "Remitting_replicate6",
"GSM4002714", "Naive_replicate5",
"GSM4002715", "Priming_replicate4",
"GSM4002716", "Priming_replicate5",
"GSM4002717", "Priming_replicate6",
"GSM4002718", "Naive_replicate6"
) %>% matrix(ncol = 2, byrow = T) %>% data.frame() %>% 
  select(GSM = X1, sample_name = X2) %>% 
  mutate(treatment = gsub( "_.*$", "",  sample_name))

dge_list <- list()
fname_vec <- list.files("0_data/other/Wheeler_GSE130119_RAW/", pattern = "*.dge.txt.gz")
for (fname in fname_vec) {
  cat(fname, "\n")
  gsm <- gsub("_.*$", "", fname)
  dge_list[[gsm]] <- read.table(gzfile(paste0("0_data/other/Wheeler_GSE130119_RAW/", fname)), header = T, row.names = 1)
}

### preprocessing and integrating
so_list <- lapply(dge_list, CreateSeuratObject, min.cells = 1) # create a list of seurat objects
saveRDS(so_list, file = paste0(dir_out, "so_list.rds"))
rm(dge_list) # free 7.4Gb memory
so_list <- lapply(names(so_list), function(i) {
  cat(i, "\n")
  x <- so_list[[i]]
  x$orig.ident <- i
  x <- subset(x, subset = nFeature_RNA > 400)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

so_anchors <- FindIntegrationAnchors(object.list = so_list, dims = 1:20)
so_combined <- IntegrateData(anchorset = so_anchors, dims = 1:20)
DefaultAssay(so_combined) <- "integrated"
saveRDS(so_combined, file = paste0(dir_out, "so_combined.rds"))

# Run the standard workflow for visualization and clustering
so_combined <- ScaleData(so_combined, verbose = FALSE)
so_combined <- RunPCA(so_combined, npcs = 15, verbose = FALSE)
# t-SNE and Clustering
so_combined <- RunUMAP(so_combined, reduction = "pca", dims = 1:15, min.dist = 0.5)
so_combined <- FindNeighbors(so_combined, reduction = "pca", dims = 1:15)
so_combined <- FindClusters(so_combined, resolution = 0.5)

#### Visualization
cluster_markers <- read.xlsx("/singerlab/linglin/Th17_single_cell_eae_ut/0_data/other/Wheeler_GSE130119_RAW/41586_2020_1999_MOESM2_ESM.xlsx", rowNames = T)
cluster_markers_ls <- lapply(0:21, function(x) {
  cluster_markers$gene[cluster_markers$cluster == x]
})
## compute and stain
DefaultAssay(so_combined) <- "RNA"
so_combined <- AddModuleScore(so_combined, features = cluster_markers_ls)

cl_ann <- c("Microglia/macrophages", 
            "Microglia/macrophages",	
            "Monocytes/macrophages",
            "Microglia/macrophages",
            "Dendritic cells",
            "Ependymal cells",
            "Oligodendrocytes",
            "T cells",
            "Astrocytes",
            "Pericytes",
            "Astrocytes",
            "Neutrophils",
            "Endothelial cells",
            "Smooth muscle cells",
            "Astrocytes",
            "Endothelial cells",
            "Neurons",
            "Epithelial cells",
            "Monocytes/macrophages",
            "Epithelial cells",
            "Stromal cells",
            "Neurons")
colnames(so_combined@meta.data)[grep("Cluster", colnames(so_combined@meta.data))] <- paste("Cluster", 0:21, cl_ann, sep = "_")

pdf(file = paste0(dir_out, "cxcl16.pdf"))
FeaturePlot(so_combined, features = c("Cxcl16"), pt.size = 0.1)
dev.off()
pdf(file = paste0(dir_out, "clusters.pdf"))
for (s in colnames(so_combined@meta.data)[grep("Cluster", colnames(so_combined@meta.data))]) {
  print(FeaturePlot(so_combined, features = s, min.cutoff = 'q5', max.cutoff = 'q95'))
}
dev.off()

saveRDS(so_combined, file = paste0(dir_out, "so_combined.rds"))
