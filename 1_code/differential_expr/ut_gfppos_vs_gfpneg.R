###################################################################################
########    Differential Expression: GFP+ vs. GFP- in each tissue (UT) ############
###################################################################################
## Identify differentially expressed genes in each tissue comparing UT GFP+ vs. UT GFP- cells
## edgeR/QLF with cellular detection rate as covariate

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

source("1_code/utils.R")
source("1_code/differential_expr/utils.R")

library(Seurat)
library(dplyr)
library(tibble)
library(edgeR)
library(openxlsx)
set.seed(1)

# configuration
today <- Sys.Date()
prep_date <- "2020-03-24"
cluster_date <- "2020-03-25"
dir_out <- paste0("2_pipeline/differential_expression/UT_GFPpos_vs_GFPneg/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_dominant_TCR_", prep_date, ".rds"))
## use clsuter annotation to filter out non-conventional Th17 cells (Treg-like and proliferating clusters)
cluster_ann <- read.xlsx(paste0("2_pipeline/clustering/UT_GFPall_intra/2020-03-25/cluster_annotation.xlsx"), sheet = "processed")
cluster_ann$tissue <- sub("_.*", "", cluster_ann$Cluster)
cluster_ann_conv <- cluster_ann %>% filter(!(Annotation %in% c("Proliferating", "Treg-like")))

########Differential Expression#########
gene_prefilter <- 0.1 # a gene should be detected by this proportion in at least one of the GFP+ or GFP- groups
### prefilter genes
prop_pos <- Matrix::rowMeans(so_all@assays$RNA@counts[,so_all$GFP_positive] > 0)
prop_neg <- Matrix::rowMeans(so_all@assays$RNA@counts[,!so_all$GFP_positive] > 0)
keep <- (prop_pos > gene_prefilter) | (prop_neg > gene_prefilter)
so_all <- subset(so_all, features = which(keep))

tissue_vec <- c("SPL", "PP", "MLN", "SI", "COL")
for (tissue in tissue_vec) {
  cat("\n=======", tissue, "=========\n")
  ### subset data
  meta_tissue <- read.table(paste0("2_pipeline/clustering/UT_GFPall_intra/", cluster_date, "/FILES/", tissue, "/meta_data.csv"), 
                            sep = ",", header = T, row.names = 1)
  conv_clusters <- as.character(cluster_ann_conv$Cluster[(cluster_ann_conv$tissue == tissue)])
  meta_tissue$tissue_cluster <- paste(meta_tissue$tissue, meta_tissue$seurat_clusters, sep = "_")
  meta_tissue <- meta_tissue[meta_tissue$tissue_cluster %in% conv_clusters,]
  so <- so_all[,rownames(meta_tissue)]
  print(dim(so))
  
  # prep data
  so$batch <- factor(so$batch)
  contrasts(so$batch) = contr.sum(nlevels(so$batch))
  so$GFP_positive <- factor(so$GFP_positive, levels = c(FALSE, TRUE))
  
  ### run edgeR
  # Create a DGEList data object.
  dgeFull <- DGEList(so[['RNA']]@counts)
  # Estimate the normalization factors
  cat("calcNormFactors\n")
  dgeFull <- calcNormFactors(dgeFull, method="TMM")
  # compute detection rate [here I use scaled number of detected genes]
  so@meta.data$cdr <- scale(so$nFeature_RNA)
  # create design matrix
  design <- model.matrix(~ GFP_positive + batch + cdr, data = so@meta.data)
  # Estimate dispersion
  cat("estimateDisp\n")
  dgeFull <- estimateDisp(dgeFull, design = design)
  # Perform QLF tests
  cat("glmQLFit&Test\n")
  fit <- glmQLFit(dgeFull, design = design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
  saveRDS(res, file = paste0(dir_out, tissue, "_results.rds"))
}
