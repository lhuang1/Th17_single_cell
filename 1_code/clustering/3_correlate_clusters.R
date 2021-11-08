####### Compute correlations for intra-tissue clusters to identify "superclusters"  ###########

# rm(list = ls())
rm(list = setdiff(ls(), "so_all"))

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(Seurat)
library(dplyr)
set.seed(1)
source("1_code/utils.R")
source("1_code/clustering/clustering_utils.R")

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
today <- Sys.Date()
clustering_date <- "2020-03-25"

treatment <- "EAE"
if (treatment == "UT") {
  tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL")
  dir_out <- paste0("2_pipeline/clustering/UT_GFPall_intra/", clustering_date, "/")
} else if (treatment == "EAE") {
  tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL", "CNS", "DLN")
  dir_out <- paste0("2_pipeline/clustering/EAE_GFPall_intra/", clustering_date, "/")
} 
stopifnot(dir.exists(dir_out))

### load data
so_all <- readRDS("2_pipeline/preprocessing/so_processed_withTCR_2020-03-16.rds")
so <- so_all[,so_all$treatment == treatment]
rm(so_all)
so@meta.data <- droplevels(so@meta.data)
DefaultAssay(so) <- "integrated"


### load clustering data
meta <- list()
markers <- list()
# tissue_cluster <- c()
for (tissue in tissue_vec) {
  meta[[tissue]]  <- read.table(paste0(dir_out, "FILES/", tissue, "/meta_data.csv"), header = T, row.names = 1, sep = ",")
  meta[[tissue]] <- meta[[tissue]][,-grep("integrated_snn_res", colnames(meta[[tissue]]))]
  markers[[tissue]] <- readRDS(paste0(dir_out, "FILES/", tissue, "/so_markers.rds"))
}

all_meta <- Reduce(rbind, meta)
so$tissue_cluster <- paste(so$tissue, all_meta[colnames(so), "seurat_clusters"], sep = "_")
tissue_cluster_vec <- unique(so$tissue_cluster)

### obtain cluster profiles
cluster_profile <- sapply(tissue_cluster_vec, function(x){
  Matrix::rowMeans(so[['integrated']]@data[so[['integrated']]@var.features, so$tissue_cluster == x])
})

### compute cosine distance
M <- as.matrix(cluster_profile) %>% t
sim <- M / sqrt(rowSums(M * M))
cosine_sim <- sim %*% t(sim)

pdf(paste0(dir_out, "cluster_similarity_", today, ".pdf"))
dist_mat <- as.dist(1 - cosine_sim)
hc <- hclust(dist_mat, method = "ward.D")
plot(hc)

hc_order = hc$order
dist_mat_transformed <- cosine_sim[hc_order, hc_order]
diag(dist_mat_transformed) <- NA

## correlation (cosine distance) plot
corrplot(dist_mat_transformed, type="upper", order="original", method = "circle",
         col=rev(brewer.pal(n=8, name="RdYlBu")),
         mar = c(5, 4, 4, 2) + 0.1,
         tl.col = "black",
         #tl.pos = 'td',
         tl.cex = 0.7,
         diag = T, na.label = " ",
         title = "Cluster Label Color by Annotation")

### dotplot
glist_df <- read.xlsx(paste0(dir_out, "cluster_annotation.xlsx"), sheet = "genes_processed")
glist <- glist_df$Gene
ct_mtx = so@assays[['RNA']]@data[glist,] %>% as.matrix()
pct_0 = sapply(colnames(dist_mat_transformed), function(x){
  rowMeans(ct_mtx[, so$tissue_cluster == x] == 0)
})

exp_mtx = so@assays[['integrated']]@data[glist,] %>% as.matrix()
ave_expr = sapply(colnames(dist_mat_transformed), function(x){
  rowMeans(exp_mtx[,so$tissue_cluster == x])
}) %>% apply(., 1, scale) %>% t

dotplot_df <- data.frame(expr = as.numeric(ave_expr), pct = 1 - as.numeric(pct_0), 
         tissue_cluster = rep(colnames(dist_mat_transformed), each = length(glist)) %>% factor(., levels = rev(colnames(dist_mat_transformed))),
         gene = rep(glist, length(colnames(dist_mat_transformed))) %>% factor(., levels = unique(glist)))
ggplot(dotplot_df) +
  geom_point(aes(x = gene, y = tissue_cluster, size = pct, color = expr)) +
  labs(x = "", y = "") +
  scale_colour_gradient(low = "lightgrey", high = "blue") +
  scale_size(range = c(0, 6), trans = scales::trans_new("square", function(x) x^2, "sqrt")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        # axis.text.y = element_text(color = rev(as.character(tissue_colors(tissue_vec)))),
        legend.title = element_text(size = 8),
        legend.box = "vertical", legend.position="right", legend.key.size = unit(.13, "in")) +
  guides(color = guide_colorbar(title = "Scaled Average\nExpression", title.vjust = .8),
         size = guide_legend(title = "Percent Cell\nDetected", title.vjust = .5, ))
dev.off()


## save data for future use
saveRDS(list("cluster_profile" = cluster_profile, 
             "cosine_sim" = cosine_sim, 
             "hierarchical_clustering" = hc, 
             "dotplot_df" = dotplot_df), 
        paste0(dir_out, "cluster_similarity_", today, ".rds"))

