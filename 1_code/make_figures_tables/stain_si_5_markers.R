#########################################################################
########   stain genes on SI GFPall UT and GFPall EAE UMAPs  ############
# rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(ggpubr)
source("1_code/utils.R")

set.seed(1)

# configuration
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-05-18"
  ut_clustering_date <- "2020-03-25"
  eae_clustering_date <- "2020-03-25"
} else {
  today <- cargs[1]
  clustering_date <- cargs[2]
}
dir_ut_clustering <- paste0("2_pipeline/clustering/UT_GFPall_intra/", ut_clustering_date, "/")
dir_eae_clustering <- paste0("2_pipeline/clustering/EAE_GFPall_intra/", eae_clustering_date, "/")
stopifnot(dir.exists(dir_ut_clustering))
stopifnot(dir.exists(dir_eae_clustering))

dir_out <- paste0("2_pipeline/SI_5/")
if (!dir.exists(dir_out)) dir.create(dir_out)



## load expression data (Seurat object)
so_all <- readRDS("2_pipeline/preprocessing/so_processed_withTCR_2020-03-16.rds")
# rm(so_all)
## load reduction and clustering (meta) data
reductions_eae <- readRDS(paste0(dir_eae_clustering, "FILES/SI/reductions.rds"))
meta_eae <- read.table(paste0(dir_eae_clustering, "FILES/SI/meta_data.csv"), sep = ",", row.names = 1, header = T)
so_eae <- so_all[,rownames(meta_eae)]
so_eae@reductions <- reductions_eae
so_eae@meta.data <- meta_eae

reductions_ut <- readRDS(paste0(dir_ut_clustering, "FILES/SI/reductions.rds"))
meta_ut <- read.table(paste0(dir_ut_clustering, "FILES/SI/meta_data.csv"), sep = ",", row.names = 1, header = T)
so_ut <- so_all[,rownames(meta_ut)]
so_ut@reductions <- reductions_ut
so_ut@meta.data <- meta_ut

reductions_both_new <- readRDS(paste0("2_pipeline/clustering/all_intra/2020-05-07/FILES/SI/reductions.rds"))
meta_both_new <- read.table(paste0("2_pipeline/clustering/all_intra/2020-05-07/FILES/SI/meta_data.csv"), sep = ",", row.names = 1, header = T)
so_new <- so_all[,rownames(meta_both_new)]
so_new@reductions <- reductions_both_new
so_new@meta.data <- meta_both_new

reductions_both_old <- readRDS(paste0("2_pipeline/clustering/all_intra/2020-04-28/FILES/SI/reductions.rds"))
meta_both_old <- read.table(paste0("2_pipeline/clustering/all_intra/2020-04-28/FILES/SI/meta_data.csv"), sep = ",", row.names = 1, header = T)
so_old <- so_all[,rownames(meta_both_old)]
so_old@reductions <- reductions_both_old
so_old@meta.data <- meta_both_old


########## stain SI_5 specific genes as a signature ########
## load DE genes
deg <- readRDS(paste0("2_pipeline/differential_expression/Cluster_EAE_GFPall_intra/2020-05-11/SI_5_results.rds"))
sig <- list("EAE_SI_5" = rownames(deg)[deg$FDR < 0.05 & deg$logFC > log2(1.5)])
so_eae <- AddModuleScore(so_eae, sig, name = "SIG")
colnames(so_eae@meta.data)[grep("SIG", colnames(so_eae@meta.data))] <- names(sig)
so_ut <- AddModuleScore(so_ut, sig, name = "SIG")
colnames(so_ut@meta.data)[grep("SIG", colnames(so_ut@meta.data))] <- names(sig)
so_old <- AddModuleScore(so_old, sig, name = "SIG")
colnames(so_old@meta.data)[grep("SIG", colnames(so_old@meta.data))] <- names(sig)
so_new <- AddModuleScore(so_new, sig, name = "SIG")
colnames(so_new@meta.data)[grep("SIG", colnames(so_new@meta.data))] <- names(sig)

p1 <- DimPlot(so_ut, group.by = "seurat_clusters", label = T)
p2 <- DimPlot(so_eae, group.by = "seurat_clusters", label = T)
q99_eae <- quantile(so_eae$EAE_SI_5, 0.99)

p3 <- data.frame(so_ut@reductions$umap@cell.embeddings) %>% 
  mutate(EAE_SI_5 = pmax(pmin(so_ut$EAE_SI_5, q99_eae), 0)) %>% 
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = EAE_SI_5), size = 0.5) +
  scale_color_gradient(low = "grey80", high = "blue", limits = c(0, ceiling(q99_eae * 10)/10)) +
  labs(color = "") +
  ggtitle("EAE_SI_5") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p4 <- data.frame(so_eae@reductions$umap@cell.embeddings) %>% 
  mutate(EAE_SI_5 = pmax(pmin(so_eae$EAE_SI_5, q99_eae), 0)) %>% 
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = EAE_SI_5), size = 0.5) +
  scale_color_gradient(low = "grey80", high = "blue", limits = c(0, ceiling(q99_eae * 10)/10)) +
  labs(color = "") +
  ggtitle("EAE_SI_5") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pdf(paste0(dir_out, "stain_SI_5_signature_", today, ".pdf"), width = 10, height = 10)
grid.arrange(p1, p2, p3, p4)
dev.off()

si_5_cells <- colnames(so_eae)[so_eae$seurat_clusters == 5]
so_old$in_EAE_SI_5 <- colnames(so_old) %in% si_5_cells
p1_old <- DimPlot(so_old, group.by = "seurat_clusters", label = T)
p2_old <- DimPlot(so_old, group.by = "in_EAE_SI_5") + scale_color_manual(values = c("grey80", "black"))
q99_old <- quantile(so_old$EAE_SI_5, 0.99)
p3_old <- data.frame(so_old@reductions$umap@cell.embeddings) %>% 
  mutate(EAE_SI_5 = pmax(pmin(so_old$EAE_SI_5, q99_old), 0),
         treatment = so_old$treatment) %>% 
  filter(treatment == "UT") %>% 
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = EAE_SI_5), size = 0.5) +
  scale_color_gradient(low = "grey80", high = "blue", limits = c(0, ceiling(q99_old * 10)/10)) +
  labs(color = "") +
  ggtitle("UT") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
p4_old <- data.frame(so_old@reductions$umap@cell.embeddings) %>% 
  mutate(EAE_SI_5 = pmax(pmin(so_old$EAE_SI_5, q99_old), 0),
         treatment = so_old$treatment) %>% 
  filter(treatment == "EAE") %>% 
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = EAE_SI_5), size = 0.5) +
  scale_color_gradient(low = "grey80", high = "blue", limits = c(0, ceiling(q99_old * 10)/10)) +
  labs(color = "") +
  ggtitle("EAE") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

so_new$in_EAE_SI_5 <- colnames(so_new) %in% si_5_cells
p1_new <- DimPlot(so_new, group.by = "seurat_clusters", label = T)
p2_new <- DimPlot(so_new, group.by = "in_EAE_SI_5") + scale_color_manual(values = c("grey80", "black"))
q99_new <- quantile(so_new$EAE_SI_5, 0.99)
p3_new<- data.frame(so_new@reductions$umap@cell.embeddings) %>% 
  mutate(EAE_SI_5 = pmax(pmin(so_new$EAE_SI_5, q99_new), 0),
         treatment = so_new$treatment) %>% 
  filter(treatment == "UT") %>% 
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = EAE_SI_5), size = 0.5) +
  scale_color_gradient(low = "grey80", high = "blue", limits = c(0, ceiling(q99_new * 10)/10)) +
  labs(color = "") +
  ggtitle("UT") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
p4_new <- data.frame(so_new@reductions$umap@cell.embeddings) %>% 
  mutate(EAE_SI_5 = pmax(pmin(so_new$EAE_SI_5, q99_new), 0),
         treatment = so_new$treatment) %>% 
  filter(treatment == "EAE") %>% 
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = EAE_SI_5), size = 0.5) +
  scale_color_gradient(low = "grey80", high = "blue", limits = c(0, ceiling(q99_new * 10)/10)) +
  labs(color = "") +
  ggtitle("EAE") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pdf(paste0(dir_out, "highlight_SI_5_", today, ".pdf"))
grid.arrange(p1_old, p2_old, p3_old, p4_old)
grid.arrange(p1_new, p2_new, p3_new, p4_new)
dev.off()


########## stain specific genes ##########
glist <- c("Cxcr6", "Itgb1", "Lgals3",  "Cd48", "Cd52")
plist_vln <- list()
plist_umap <- list()
plist_cluster <- list()
for (treatment in c("EAE", "UT")) {
  if (treatment == "EAE") {
    so <- so_eae
  } else {
    so <- so_ut
  }
  DefaultAssay(so) <- "integrated"
  plist_cluster[[treatment]] <- DimPlot(so, group.by = "seurat_clusters", label = T)
  
  plist_vln[[treatment]] <- lapply(glist, function(g) {
    cat(g, "\n")
    if (g %in% rownames(so[["integrated"]])) {
      p <- VlnPlot(so, features = g, group.by = "seurat_clusters")
    } else {
      p <- ggplot()
    }
    p
  })
  names(plist_vln[[treatment]]) <- glist
  plist_umap[[treatment]] <- lapply(glist, function(g) {
    cat(g, "\n")
    if (g %in% rownames(so[["integrated"]])) {
      # p <- FeaturePlot(so, features = g, min.cutoff = "q1", max.cutoff = "q99")
      q99_both <- quantile(c(as.numeric(so_ut[["integrated"]][g,], so_eae[["integrated"]][g,])), 0.99)
      
      p <- data.frame(so@reductions$umap@cell.embeddings) %>% 
        mutate(expr = pmax(pmin(as.numeric(so[["integrated"]][g,]), q99_both), 0)) %>% 
        ggplot() +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = expr), size = 0.2) +
        scale_color_gradient(low = "grey80", high = "blue", limits = c(0, ceiling(q99_both * 10)/10)) +
        labs(color = "") +
        ggtitle(g) +
        theme_classic() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    } else {
      p <- ggplot()
    }
    p
  })
  names(plist_umap[[treatment]]) <- glist
  
}

pdf(paste0(dir_out, "stain_SI_5_markers_", today, ".pdf"))
print(plist_cluster$EAE)
print(plist_cluster$EAE)
for (g in glist) {
  cat(g, "\n")
  plist <- list(plist_vln$UT[[g]], plist_vln$EAE[[g]], plist_umap$UT[[g]], plist_umap$EAE[[g]])
  grid.arrange(grobs = plist)
}
dev.off()



### bulk signature
sig_bulk <- readRDS("2_pipeline/SI_5/bulk/2020-05-11/signature_cd29_pos_vs_neg_filtered.rds")
sapply(sig_bulk, length)
sig_bulk_intersect <- sapply(sig_bulk, function(x){intersect(x, rownames(so_eae))})
sapply(sig_bulk_intersect, length)
so_eae <- AddModuleScore(so_eae, features = sig_bulk_intersect, name = "BULK")
colnames(so_eae@meta.data)[grep("BULK", colnames(so_eae@meta.data))] <- names(sig_bulk_intersect)

pdf(paste0(dir_out, "stain_bulk_", today, ".pdf"))
for (g in names(sig_bulk_intersect)) {
  p <- VlnPlot(so_eae, features = g, group.by = "seurat_clusters")
  print(p)
  p <- FeaturePlot(so_eae, features = g, min.cutoff = "q1", max.cutoff = "q99")
  print(p)
}
dev.off()
