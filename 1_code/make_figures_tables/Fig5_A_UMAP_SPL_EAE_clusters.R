#### visualize 5 SPL EAE clusters on UMAP
rm(list = setdiff(ls(), "so_all"))

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
source("1_code/utils.R")

#### configuration ####
today <- "2020-11-02"
clustering_date <- "2020-03-25"

dir_out <- paste0("3_output/Figure_5/")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

#### load data ####
meta <- read.table(paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/FILES/meta_data_switched.txt"), 
                   header = T, sep = " ")
reductions <- readRDS(paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/FILES/reductions.rds"))

#### make plot ####
stopifnot(all(rownames(meta) == rownames(reductions$umap@cell.embeddings)))
plt_df <- data.frame(reductions$umap@cell.embeddings)
plt_df$cluster <- meta$seurat_clusters

p <- ggplot(plt_df, aes(x = UMAP_1, y = UMAP_2, color = factor(cluster))) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = SPL_EAE_all_cluster_colors()) +
  labs(color = "Cluster") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))

pdf(file = paste0(dir_out, "Fig5_A_UMAP_SPL_EAE_clusters.pdf"), width = 7, height = 7)
p
dev.off()