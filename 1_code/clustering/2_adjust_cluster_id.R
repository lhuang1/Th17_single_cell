##### Switch cluster IDs for SPL_EAE and CNS_EAE to make numbering meaningful ######
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

## Switch EAE SPL_0 and SPL_1, such that SPL_0 is non-pathogenic (Slamf6+) and SPL_1 is pathogenic (Cxcr6+)
file.rename(from = "2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/SPL/meta_data.csv",
            to = "2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/SPL/meta_data_before_switch.csv")
meta <- read.table("2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/SPL/meta_data_before_switch.csv", sep = ",", header = T)
new_cluster <- meta$seurat_clusters
new_cluster[new_cluster %in% 0:1] <- 1 - new_cluster[new_cluster %in% 0:1]
table(new_cluster, meta$seurat_clusters)
meta$seurat_clusters <- new_cluster
write.table(meta, paste0("2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/SPL/meta_data.csv"), sep = ",")


## Switch UT SPL_0 and SPL_1, such that SPL_0 is non-pathogenic (Slamf6+) 
file.rename(from = "2_pipeline/clustering/UT_GFPall_intra/2020-03-25/FILES/SPL/meta_data.csv",
            to = "2_pipeline/clustering/UT_GFPall_intra/2020-03-25/FILES/SPL/meta_data_before_switch.csv")
meta <- read.table("2_pipeline/clustering/UT_GFPall_intra/2020-03-25/FILES/SPL/meta_data_before_switch.csv", sep = ",", header = T, row.names = 1)
new_cluster <- meta$seurat_clusters
new_cluster[new_cluster %in% 0:1] <- 1 - new_cluster[new_cluster %in% 0:1]
table(new_cluster, meta$seurat_clusters)
meta$seurat_clusters <- new_cluster
write.table(meta, paste0("2_pipeline/clustering/UT_GFPall_intra/2020-03-25/FILES/SPL/meta_data.csv"), sep = ",")

## Switch CNS clusters to match pseudo-time ordering
file.rename(from = "2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/CNS/meta_data.csv",
            to = "2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/CNS/meta_data_before_switch.csv")
meta <- read.table("2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/CNS/meta_data_before_switch.csv", header = T, sep = ",", row.names = 1)
### switch cluster labels!!!
# Quiescent: C1 -> C0
# Migratory: C0 -> C1
# Proliferating: C3 -> C2
# CD8-like: C4 -> C3
# Proinflammatory: C2 -> C4
cluster_new <- meta$seurat_clusters
cluster_new[meta$seurat_clusters == 1] <- 0
cluster_new[meta$seurat_clusters == 0] <- 1
cluster_new[meta$seurat_clusters == 3] <- 2
cluster_new[meta$seurat_clusters == 4] <- 3
cluster_new[meta$seurat_clusters == 2] <- 4
meta$seurat_clusters <- cluster_new
write.table(meta, "2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/CNS/meta_data.csv", sep = ",")
