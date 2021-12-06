######## cluster cells across multiple tissues ############

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
set.seed(1)
source("1_code/utils.R")
source("1_code/clustering/utils.R")

#### configuration ####
param_ls <- list(treatment = "UT", gfp = "GFPpos", PC_dims = 1:20, cluster_res = 0.3, UMAP_nneighbors = 30, UMAP_mindist = 0.3)

today <- Sys.Date()
dir_out <- paste0("2_pipeline/clustering/", param_ls$treatment, "_", param_ls$gfp, "_inter/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
  dir.create(paste0(dir_out, "PLOTS"))
  dir.create(paste0(dir_out, "FILES"))
}


### load data
prep_type_date <- "dominant_TCR_2020-03-24"
so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_", prep_type_date, ".rds"))
so <- so_all[,so_all$treatment == param_ls$treatment]
if (param_ls$gfp == "GFPpos") {
  so <- so[,so$GFP_positive]
} else if (param_ls$gfp == "GFPneg") {
  so <- so[,!so$GFP_positive]
} else if (param_ls$gfp != "GFPall"){
  stop("Double check param_ls$gfp.")
}
rm(so_all)
DefaultAssay(object = so) <- "integrated"


### clustering
# PCA
so <- RunPCA(object = so)
pdf(paste0(dir_out, "/PLOTS/PCA_Elbow.pdf"))
print(ElbowPlot(object = so, ndims = 50))
dev.off()
write.csv(Embeddings(so, reduction = "pca"), file = paste0(dir_out, "/FILES/coords_pca.csv"), quote = F)

# clustering
so <- FindNeighbors(object = so, dims = param_ls$PC_dims)
so <- FindClusters(object = so, resolution = param_ls$cluster_res)
write.csv(so@meta.data, file = paste0(dir_out, "/FILES/meta_data.csv"))

# UMAP
so <- RunUMAP(object = so, verbose = TRUE, 
              dims = param_ls$PC_dims, n.neighbors = param_ls$UMAP_nneighbors, min.dist = param_ls$UMAP_mindist) 
write.csv(Embeddings(so, reduction = "umap"), file = paste0(dir_out, "/FILES/coords_umap.csv"), quote = F)

# overview and cluster compositions (pie chart)
pt.size <- 0.1
pie_by <- c("tissue", "GFP_positive", "batch")
cols_by <- c("seurat_clusters", "tissue", "batch", "orig.ident", "GFP_positive")
plot_cluster_overview(so, dir_out, pt.size = pt.size, pie_by = pie_by, cols_by = cols_by)

### staining ###
### load genes of interest
genes_of_interest <- read.table("0_data/gene_lists/genes_of_interest_Alex_20190421.csv", stringsAsFactors = FALSE)[,1]
genes_of_interest <- c(genes_of_interest, "Il17a-GFP", "TdTomato")

# genes of interest
pdf(file = paste0(dir_out, "/PLOTS/UMAP_genes_of_interest.pdf"))
for (g in genes_of_interest) {
  print(g)
  if (! (g %in% rownames(so))) {
    cat(g, ": not in data...\n"); next
  }
  print(VlnPlot(object = so, features = g, pt.size = pt.size))
  print(FeaturePlot(object = so, features = g,
                    min.cutoff = "q5", max.cutoff = "q95", reduction = "umap", pt.size = pt.size))
}
dev.off()

# ### Identify and plot cluster marker genes ###
so_markers <- FindAllMarkers(object = so, only.pos = TRUE, min.pct = 0.3, max.cells.per.ident = 500,
                             logfc.threshold = 0.25)
saveRDS(so_markers, paste0(dir_out, "/FILES/so_markers.rds"))
save_and_plot_markers(so, so_markers, dir_markers = dir_out, pt.size = pt.size)


