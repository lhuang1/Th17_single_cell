#### Clean out non-Th17 cells and batch-specific clusters #####

#### configuration ####
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-01-29"
  integrate_date <- "2020-01-29"
  clean_round <- "3"
  run_clustering <- TRUE
} else {
  today <- cargs[1]
  integrate_date <- cargs[2]
  clean_round <- cargs[3]
  run_clustering <- as.logical(cargs[4])
}
dir_out <- paste0("2_pipeline/preprocessing/cleaning/", today, "/round_", clean_round, "/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)


library(Seurat)
library(ggplot2)
library(dplyr)
source("1_code/utils.R")
source("1_code/preprocessing/utils.R")

so_integrated_all <- readRDS(file = paste0("2_pipeline/preprocessing/integration/", integrate_date, "/so_integrated_SCTransformed.rds"))

if (clean_round != "0") { # not initial round; need to take out clusters
  # read clusters to take out
  dir_last_round <- paste0("2_pipeline/preprocessing/cleaning/", today, "/round_", as.numeric(clean_round) - 1, "/")
  take_out <- read.csv(file = paste0(dir_last_round, "clusters_take_out.csv"), header = F, row.names = 1)
  keep_cells <- list()
  for (i in rownames(take_out)) {
    clusters_i <- as.numeric(take_out[i,])
    clusters_i <- clusters_i[!is.na(clusters_i)]
    meta <- read.table(file = paste0(dir_last_round, i, "/FILES/meta_data.txt"), stringsAsFactors = F)
    keep_cells[[i]] <- rownames(meta)[!meta$seurat_clusters %in% clusters_i]
  }
  so_integrated_all <- so_integrated_all[,unlist(keep_cells)]
  saveRDS(colnames(so_integrated_all), file = paste0(dir_out, "keep_cell_names.rds"))
} 

tissue_vec <- c("SPL", "PP", "MLN", "SI", "COL", "CNS", "DLN")
treatment_vec <- c("UT", "EAE")
# load genes of interest and contamination markers
cont_markers <- read.table("0_data/gene_lists/Genes_for_cleanup.csv", stringsAsFactors = FALSE)[,1]
cont_markers <- c(cont_markers, "Ptprc", "Tcrg-C4", "Tcrg-C2", "Tcrg-C1", "Il17a-GFP") ## add a few genes
genes_of_interest <- read.table("0_data/gene_lists/genes_of_interest_Alex_20190421.csv", stringsAsFactors = FALSE)[,1]
genes_of_interest <- c(genes_of_interest, "Il17a-GFP")

batch_clusters <- list()

for (tissue in tissue_vec) {
  for (treatment in treatment_vec) {
    if (treatment == 'UT' & tissue %in% c("CNS", "DLN")) next
    sample_name <- paste(tissue, treatment, sep = "_")
    cat("\n\n", sample_name, "...\n")
    dir.create(paste0(dir_out, sample_name, "/FILES"), recursive = T, showWarnings = F)
    dir.create(paste0(dir_out, sample_name, "/PLOTS"), recursive = T, showWarnings = F)
    
    so_integrated <- so_integrated_all[,so_integrated_all$tissue == tissue & so_integrated_all$treatment == treatment]
    
    ## Clustering
    so_integrated <- ScaleData(object = so_integrated, verbose = TRUE)
    so_integrated <- RunPCA(object = so_integrated, npcs = 50, verbose = TRUE)
    pdf(paste0(dir_out, sample_name, "/PLOTS/PCA_Elbow.pdf"))
    print(ElbowPlot(object = so_integrated, ndims = 50))
    dev.off()
    write.table(Embeddings(so_integrated, reduction = "pca"), 
                file = paste0(dir_out, sample_name, "/FILES/coords_pca.txt"), quote = F)
    
    # clusteirng
    so_integrated <- FindNeighbors(object = so_integrated, dims = 1:20)
    so_integrated <- FindClusters(object = so_integrated, resolution = 1.2)
    write.table(so_integrated@meta.data, 
                file = paste0(dir_out, sample_name, "/FILES/meta_data.txt"), quote = F)
    
    # UMAP
    so_integrated <- RunUMAP(object = so_integrated, verbose = FALSE, dims = 1:20, 
                             n.neighbors = 30, min.dist = 0.3, seed.use = 1)
    write.table(Embeddings(so_integrated, reduction = "umap"), 
                file = paste0(dir_out, sample_name, "/FILES/coords_umap.txt"), quote = F)
    
    # plot overview
    pdf(file = paste0(dir_out, sample_name, "/PLOTS/UMAP_overview.pdf"), width = 12)
    p <- so_integrated@meta.data %>%
      group_by(seurat_clusters) %>% summarize(count = n()) %>%
      ggplot(aes(x=seurat_clusters, y = count))+
      geom_bar(aes(fill = seurat_clusters), stat = "identity") +
      geom_text(aes(label=count), vjust=0) +
      ggtitle("Cluster sizes")
    print(p)
    pt.size = 0.2
    print(DimPlot(object = so_integrated, reduction = "umap", group.by = "seurat_clusters",
                  pt.size = pt.size,
                  label = TRUE))
    print(DimPlot(object = so_integrated, reduction = "umap", group.by = "batch",
                  pt.size = pt.size))
    print(DimPlot(object = so_integrated, reduction = "umap", group.by = "tissue",
                  pt.size = pt.size))
    print(DimPlot(object = so_integrated, reduction = "umap", group.by = "treatment", 
                  pt.size = pt.size))
    print(DimPlot(object = so_integrated, reduction = "umap", group.by = "GFP_positive", 
                  cols = c("lightgrey", "green"), pt.size = pt.size))
    p <- so_integrated@meta.data %>%
      group_by(orig.ident) %>% summarize(count = n()) %>%
      ggplot(aes(x=orig.ident, y = count))+
      geom_bar(aes(fill = orig.ident), stat = "identity") +
      geom_text(aes(label=count), vjust=0) +
      ggtitle("Sample sizes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
    dev.off()
    
    DefaultAssay(object = so_integrated) <- "integrated"
    # contamination markers
    pdf(file = paste0(dir_out, sample_name, "/PLOTS/UMAP_contamination.pdf"), width = 12)
    for (g in cont_markers) {
      print(g)
      if (! (g %in% rownames(so_integrated))) {
        cat(g, ": not in data...\n"); next
      }
      print(VlnPlot(object = so_integrated, features = g, pt.size = pt.size))
      print(FeaturePlot(object = so_integrated, features = g,
                        min.cutoff = "q5", reduction = "umap", pt.size = pt.size))
    }
    dev.off()
    
    # genes of interest
    pdf(file = paste0(dir_out, sample_name, "/PLOTS/UMAP_genes_of_interest.pdf"), width = 12)
    for (g in genes_of_interest) {
      print(g)
      if (! (g %in% rownames(so_integrated))) {
        cat(g, ": not in data...\n"); next
      }
      print(VlnPlot(object = so_integrated, features = g, pt.size = pt.size))
      print(FeaturePlot(object = so_integrated, features = g,
                        min.cutoff = "q5", reduction = "umap", pt.size = pt.size))
    }
    dev.off()
    
    ## Identify markers for each cluster
    so_markers <- FindAllMarkers(object = so_integrated, 
                                 only.pos = TRUE, min.pct = 0.25, max.cells.per.ident = 500,
                                 logfc.threshold = 0.25)
    saveRDS(so_markers, paste0(dir_out, sample_name, "/FILES/so_markers.rds"))
    save_and_plot_markers(so_integrated, so_markers, dir_markers = paste0(dir_out, sample_name))
    
    ## Identify batch specific clusters
    batch_clusters[[sample_name]] <- find_batch_cluster(so_integrated@meta.data)
  }
}

saveRDS(batch_clusters, file = paste0(dir_out, "batch_specific_clusters.rds"))
