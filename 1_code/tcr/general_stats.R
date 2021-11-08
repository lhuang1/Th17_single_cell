### 

today <- Sys.Date()#"2020-01-19"
dir_proj <- "/singerlab/linglin/Th17_single_cell_eae_ut/"
dir_out <- "2_pipeline/TCR/"

library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(openxlsx)

tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL", "CNS", "DLN")
treatment_vec <- c("UT", "EAE")
### load clustering data
meta <- list()
reductions <- list()
clustering_date <- "2020-03-25"
for (treatment in treatment_vec) {
  for (tissue in tissue_vec) {
    tt <- paste(treatment, tissue, sep = "_")
    if (tt %in% c("UT_CNS", "UT_DLN")) next
    if (tt == "EAE_SPL") {
      dir_spl_cluster <- paste0(dir_proj, "2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/")
      meta[[tt]]  <- read.table(paste0(dir_spl_cluster, "FILES/meta_data_switched.txt"), header = T, sep = " ")
      meta[[tt]] <- meta[[tt]][,-grep("integrated_snn_res", colnames(meta[[tt]]))]
      meta[[tt]]$tissue_cluster <- paste(meta[[tt]]$tissue, meta[[tt]]$seurat_clusters, sep = "_")
      reductions[[tt]] <- readRDS(paste0(dir_spl_cluster, "FILES/reductions.rds"))
    } else {
      dir_clustering <- paste0(dir_proj, "2_pipeline/clustering/", treatment, "_GFPall_intra/", clustering_date, "/")
      meta[[tt]]  <- read.table(paste0(dir_clustering, "FILES/", tissue, "/meta_data.csv"), header = T, row.names = 1, sep = ",")
      meta[[tt]] <- meta[[tt]][,-grep("integrated_snn_res", colnames(meta[[tt]]))]
      meta[[tt]]$tissue_cluster <- paste(meta[[tt]]$tissue, meta[[tt]]$seurat_clusters, sep = "_")
      reductions[[tt]] <- readRDS(paste0(dir_clustering, "FILES/", tissue, "/reductions.rds"))
    }
    if (tt == "EAE_MLN") { ## take out contaminating cells ad-hoc; need to fix up-stream
      idx_contminated <- which(reductions[[tt]]$umap@cell.embeddings[,"UMAP_2"] > 8)
      reductions[[tt]]$umap@cell.embeddings <- reductions[[tt]]$umap@cell.embeddings[-idx_contminated,]
      meta[[tt]] <- meta[[tt]][-idx_contminated,]
    }
  }
}
all_meta <- Reduce(rbind, meta)
all_meta$tissue_cluster <- paste(all_meta$tissue, all_meta$seurat_clusters, sep = "_")


clonotyping_type_date <- "dominant_2020-04-02"
clty_all <- read.csv(paste0(dir_proj, "2_pipeline/TCR/clonotype_assignment_", clonotyping_type_date, ".csv"), 
                     row.names = 1)

## merge clonotypes IDs to meta data
all_meta$clonotype_id <- clty_all[rownames(all_meta),"clonotype_id"]
meta_noNA <- all_meta[!is.na(all_meta$clonotype_id),] %>% droplevels()
meta_noNA$pooled_cid <- paste(meta_noNA$treatment, meta_noNA$batch, meta_noNA$mouse, meta_noNA$clonotype_id, sep = "_")

#### tabulate number of cells/clonotypes by tissue (row) and treatment/mouse ####
df_cell_ls <- list()
df_clone_ls <- list()
for (trt in treatment_vec) {
  meta_noNA_sub <- meta_noNA %>% filter(treatment == trt)
  meta_noNA_sub$sample <- paste(meta_noNA_sub$treatment, meta_noNA_sub$batch, meta_noNA_sub$mouse, sep = "_") %>% 
    factor(labels = paste(trt, 1:4, sep = "_"))
  ## cell number
  df_cell <- with(meta_noNA_sub, table(tissue, sample)) %>%
    data.frame() %>% spread(key = sample, value = Freq) %>% 
    column_to_rownames(var = "tissue")
  df_cell <- rbind(df_cell, total = colSums(df_cell))
  df_cell <- cbind(df_cell, total = rowSums(df_cell))
  df_cell_ls[[trt]] <- df_cell
  
  ## clonotype number
  df_clone <- with(meta_noNA_sub[!duplicated(paste0(meta_noNA_sub$pooled_cid, meta_noNA_sub$tissue)),], 
                   table(tissue, sample)) %>%
    data.frame() %>% spread(key = sample, value = Freq) %>% 
    column_to_rownames(var = "tissue")
  df_clone <- rbind(df_clone, total = with(meta_noNA_sub[!duplicated(meta_noNA_sub$pooled_cid),], table(sample)) %>% as.numeric())
  
  meta_noNA_sub$TRA_TRB <- paste(meta_noNA_sub$TRA, meta_noNA_sub$TRB, sep = "_")
  df_clone <- cbind(df_clone, total = c(with(meta_noNA_sub[!duplicated(meta_noNA_sub$TRA_TRB),], 
                                             table(tissue)) %>% as.numeric(),
                                        length(unique(meta_noNA_sub$TRA_TRB))))
  df_clone_ls[[trt]] <- df_clone
}

df_cell_ut_eae <- cbind(df_cell_ls$UT, df_cell_ls$EAE)
df_clone_ut_eae <- cbind(df_clone_ls$UT, df_clone_ls$EAE)

write.csv(df_cell_ut_eae, quote = F, file = paste0(dir_out, "n_cell_by_tissue_sample", today, ".csv"))
write.csv(df_clone_ut_eae, quote = F, file = paste0(dir_out, "n_clone_by_tissue_sample", today, ".csv"))



