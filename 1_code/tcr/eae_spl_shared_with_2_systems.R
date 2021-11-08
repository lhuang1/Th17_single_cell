today <- Sys.Date()#"2020-05-04"
clustering_date <- "2020-03-25"
clonotyping_type_date <- "dominant_2020-04-02"
prep_type_date <- "dominant_TCR_2020-03-24"
trt_gfp <- "EAE_GFPall"
treatment <- gsub("_.*", "", trt_gfp)
gfp <- gsub(".*_", "", trt_gfp)
tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL", "CNS", "DLN")
dir_cluster <- paste0("2_pipeline/clustering/EAE_GFPall_intra/", clustering_date, "/")
dir_spl_cluster <- paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/")

meta <- list()
reductions <- list()
for (tissue in tissue_vec) {
  if (tissue == "SPL" & treatment == "EAE") {
    meta[[tissue]]  <- read.table(paste0(dir_spl_cluster, "FILES/meta_data_switched.txt"), header = T, sep = " ")
    meta[[tissue]] <- meta[[tissue]][,-grep("integrated_snn_res", colnames(meta[[tissue]]))]
    reductions[[tissue]] <- readRDS(paste0(dir_spl_cluster, "FILES/reductions.rds"))
  } else {
    meta[[tissue]]  <- read.table(paste0(dir_cluster, "FILES/", tissue, "/meta_data.csv"), header = T, row.names = 1, sep = ",")
    meta[[tissue]] <- meta[[tissue]][,-grep("integrated_snn_res", colnames(meta[[tissue]]))]
    reductions[[tissue]] <- readRDS(paste0(dir_cluster, "FILES/", tissue, "/reductions.rds"))
  }
}
meta_all <- Reduce(rbind, meta)
##$$$ fix batches 
meta_all$batch <- as.character(meta_all$batch)
meta_all$batch[meta_all$treatment == "UT" & meta_all$batch == "b8"] <- "b9"
meta_all$batch[meta_all$treatment == "UT" & meta_all$batch == "b7"] <- "b8"
meta_all$batch <- factor(meta_all$batch)


## clonotypes
clty_all <- read.csv(paste0("2_pipeline/TCR/clonotype_assignment_", clonotyping_type_date, ".csv"), row.names = 1)

## cluster annotations
cluster_ann <- read.xlsx(paste0(dir_cluster, "/cluster_annotation.xlsx"), sheet = "processed")
cluster_ann$Annotation[is.na(cluster_ann$Annotation)] <- "Other"
cluster_ann$Annotation <- factor(cluster_ann$Annotation, levels = unique(cluster_ann$Annotation))

## merge clonotype info to meta data
meta_all$clonotype_id <- clty_all[rownames(meta_all),"clonotype_id"]
meta_all$pooled_clonotype_id <- paste(meta_all$treatment, meta_all$batch, meta_all$mouse, meta_all$clonotype_id, sep = "_")
meta_all$tissue_cluster <- paste(meta_all$tissue, meta_all$seurat_clusters, sep = "_")
meta_all$sample <- paste(meta_all$batch, meta_all$mouse, sep = "_") %>% factor()
meta_noNA <- meta_all[!is.na(meta_all$clonotype_id),] %>% droplevels()

stopifnot(all(meta_noNA$treatment == treatment)) ## check if treatment is correct (should be UT cells only)

tissue_1 <- c("CNS", "DLN")
tissue_2 <- c("MLN", "PP", "SI", "COL")
df <- meta_noNA %>% group_by(pooled_clonotype_id) %>%
  mutate(in_CNS_DLN = any(tissue %in% tissue_1), 
         in_MLN_PP_SI_COL = any(tissue %in% tissue_2)) %>% 
  ungroup() %>% 
  filter(tissue == "SPL")
table(df$in_CNS_DLN, df$in_MLN_PP_SI_COL) / nrow(df)

table(df$in_CNS_DLN) / nrow(df)

table(df$in_CNS_DLN, df$in_MLN_PP_SI_COL) / nrow(df)
table(df$in_CNS_DLN, df$in_MLN_PP_SI_COL) / sum(df$in_CNS_DLN | df$in_MLN_PP_SI_COL)
