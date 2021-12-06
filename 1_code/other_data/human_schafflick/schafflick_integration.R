####### Integrate human (Schafflick) data ########
# 42,969 blood single-cell transcriptomes (five control vs. five MS donors) 
# and 22,357 corresponding CSF single-cell transcriptomes (four control vs. four MS donors).


### configuration
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(openxlsx)
source("1_code/utils.R")

options(future.globals.maxSize = +Inf)


today <- "2020-06-01"
dir_out <- paste0("2_pipeline/human_shafflick/", today, "/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

dir_raw <- "0_data/other/Schafflick_GSE130119_RAW/"
fnames_vec <- list.files(dir_raw, pattern = "*.gz")

samples_csf <- c('MS19270_CSF', 'MS49131_CSF', 'MS58637_CSF', 'MS71658_CSF', 'PTC32190_CSF', 'PTC41540_CSF', 'PST45044_CSF', 'PTC85037_CSF')
samples_pbmcs <- c('MS19270_PBMCs', 'MS49131_PBMCs','MS71658_PBMCs','MS60249_PBMCs','MS74594_PBMCs', 'PTC32190_PBMCs','PTC41540_PBMCs','PTC85037_PBMCs','PST83775_PBMCs','PST95809_PBMCs')

samples_all <- c(samples_csf, samples_pbmcs)

## load raw data
data_raw <- list()
for (data_name in samples_all) {
  cat("Reading", data_name, "...\n")
  gsm <- gsub("_.*$", "", fnames_vec[grep(data_name, fnames_vec)[1]])
  mtx <- Matrix::readMM(file = paste0(dir_raw, gsm, "_", data_name, "_GRCh38_matrix.mtx.gz"))
  cell_names <- readLines(paste0(dir_raw, gsm, "_", data_name, "_GRCh38_barcodes.tsv.gz"))
  colnames(x = mtx) <- paste0(data_name, "_", cell_names)
  feature_names <- read.table(gzfile(paste0(dir_raw, gsm, "_", data_name, "_GRCh38_genes.tsv.gz")), stringsAsFactors = F)[,2]
  rownames(x = mtx) <- feature_names
  data_raw[[data_name]] <- mtx
}

## convert to Seurat object and normalize
for (i in names(data_raw)) {
  cat("\n\n", i, "\n")
  x <- CreateSeuratObject(data_raw[[i]])
  data_raw[[i]] <- NULL
  x$orig.ident <- i
  x <- SCTransform(x)
  saveRDS(x, paste0(dir_out, i, "_sctransformed.rds"))
}

## integrate data
so_list <- lapply(samples_all, function(i){
  readRDS(paste0(dir_out, i, "_sctransformed.rds"))
})
names(so_list) <- samples_all
so_features <- SelectIntegrationFeatures(object.list = so_list, nfeatures = Inf)
so_list <- PrepSCTIntegration(object.list = so_list, anchor.features = so_features,
                              verbose = FALSE)
so_anchors_CSF <- FindIntegrationAnchors(object.list = so_list[grep("CSF", names(so_list))], normalization.method = "SCT",
                                     anchor.features = so_features, verbose = TRUE)
so_anchors_PBMCs <- FindIntegrationAnchors(object.list = so_list[grep("PBMCs", names(so_list))], normalization.method = "SCT",
                                         anchor.features = so_features, verbose = TRUE)
rm(so_list)

so_combined_CSF <- IntegrateData(anchorset = so_anchors_CSF, normalization.method = "SCT", verbose = TRUE)
rm(so_anchors_CSF)
saveRDS(so_combined_CSF, file = paste0(dir_out, "so_combined_CSF.rds"))

so_combined_PBMCs <- IntegrateData(anchorset = so_anchors_PBMCs, normalization.method = "SCT",verbose = TRUE)
rm(so_anchors_PBMCs)
saveRDS(so_combined_PBMCs, file = paste0(dir_out, "so_combined_PBMCs.rds"))


##### Clustering #########
## load DE genes given in the paper
f_de <- "0_data/other/Schafflick_GSE130119_RAW/41467_2019_14118_MOESM2_ESM.xlsx"
cluster_name_vec <- getSheetNames(f_de)
cluster_markers <- lapply(cluster_name_vec, function(cluster_name) {
  tmp <- read.xlsx("0_data/other/Schafflick_GSE130119_RAW/41467_2019_14118_MOESM2_ESM.xlsx", rowNames = T, sheet = cluster_name)
  tmp %>% tibble::rownames_to_column(var = "gene") %>%
    arrange(desc(logFC)) %>% filter(fdr_edgeR < 0.05) %>%
    select(gene) %>% unlist %>% as.character()
})
names(cluster_markers) <- cluster_name_vec
n <- 200 ## only use top 200 up-regulated genes
cluster_markers_n <- lapply(cluster_markers, function(x){x[1:min(n, length(x))]})

## load signatures
c7_sig_human <- readRDS("2_pipeline/other/c7_human_signatures.rds")
th17_sig_names_human <- c("GSE32901_NAIVE_VS_TH17_ENRICHED_CD4_TCELL_DN", "GSE32901_TH17_EMRICHED_VS_TH17_NEG_CD4_TCELL_UP", "GSE32901_TH1_VS_TH17_ENRICHED_CD4_TCELL_DN")
th17_sig <- c7_sig_human[th17_sig_names_human]


tissue_vec <- c("PBMCs", "CSF")
for (tissue in tissue_vec) {
  so_combined <- readRDS(file = paste0(dir_out, "so_combined_", tissue, ".rds"))
  DefaultAssay(so_combined) <- "integrated"
  
  so_combined$individual <- gsub("_.*$", "", so_combined$orig.ident)
  so_combined$tissue <- tissue #gsub(".*_", "", so_combined$orig.ident)
  so_combined$is_MS <- ifelse(grepl("MS", so_combined$individual), "MS", "Healthy")

  # Run the standard workflow for visualization and clustering
  so_combined <- RunPCA(so_combined, npcs = 15, verbose = FALSE)
  # t-SNE and Clustering
  so_combined <- RunUMAP(so_combined, reduction = "pca", dims = 1:15, min.dist = 0.5)
  so_combined <- FindNeighbors(so_combined, reduction = "pca", dims = 1:15)
  so_combined <- FindClusters(so_combined, resolution = 0.5)

  # compute DE gene score and stain on umap
  so_combined <- AddModuleScore(so_combined, features = cluster_markers_n, search = F, assay = "RNA", nbin = 18)
  colnames(so_combined@meta.data)[grep("Cluster", colnames(so_combined@meta.data))] <- cluster_name_vec
  
  so_combined <- AddModuleScore(so_combined, features = th17_sig, search = F, assay = "RNA", nbin = 18)
  colnames(so_combined@meta.data)[grep("Cluster", colnames(so_combined@meta.data))] <- names(th17_sig)
  
  pdf(file = paste0(dir_out, "cell_type_top200logFC_", tissue, ".pdf"))
  DimPlot(so_combined, label = T)
  for (s in cluster_name_vec) {
    print(FeaturePlot(so_combined, features = s, min.cutoff = 'q5', max.cutoff = 'q95'))
  }
  print(FeaturePlot(so_combined, features = "RORA", min.cutoff = 'q5', max.cutoff = 'q95'))
  dev.off()
  
  pdf(file = paste0(dir_out, "th17_signatures_", tissue, ".pdf"))
  DimPlot(so_combined, label = T)
  for (s in names(th17_sig)) {
    print(FeaturePlot(so_combined, features = s, min.cutoff = 'q5', max.cutoff = 'q95'))
  }
  dev.off()
  saveRDS(so_combined@meta.data, file = paste0(dir_out, "meta_combined_clustered_", tissue, ".rds"))
  saveRDS(so_combined@reductions, file = paste0(dir_out, "reductions_combined_clustered_", tissue, ".rds"))
}

  