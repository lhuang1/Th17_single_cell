so <- readRDS("2_pipeline/cxcl16/2020-02-26/so_combined.rds")
DimPlot(so, label = T)

so_T <- so[,so$seurat_clusters == 10]

so_T@meta.data <- droplevels(so_T@meta.data)

DefaultAssay(so_T) <- "integrated"
so_T <- RunPCA(so_T)
so_T <- RunUMAP(so_T, dims = 1:10, n.neighbors = 10, min.dist = 0.01)
so_T <- FindNeighbors(so_T)
so_T <- FindClusters(so_T, resolution = 0.3)

DimPlot(so_T, label = T)
FeaturePlot(so_T, features = "Cd8a")
VlnPlot(so_T, features = "Cd8a")
FeaturePlot(so_T, features = "Cd4")
VlnPlot(so_T, features = "Cd4")


rownames(so_T[['RNA']])[grep("Cd8", rownames(so_T[['RNA']]))]


