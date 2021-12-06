######################################################################################
########     Differential Expression: EAE vs. UT (Stratify by Tissue)     ############
######################################################################################
## Identify differentially expressed genes in each tissue comparing EAE vs. UT GFP all
## edgeR/QLF with cellular detection rate as covariate

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

source("1_code/utils.R")
source("1_code/differential_expr/utils.R")

library(Seurat)
library(dplyr)
library(edgeR)
set.seed(1)

# configuration
today <- Sys.Date()
prep_date <- "2020-03-24"
dir_out <- paste0("2_pipeline/differential_expression/EAE_vs_UT_GFPall/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_dominant_TCR_", prep_date, ".rds"))

########Differential Expression#########
tissue_vec <- c("SPL", "COL", "SI", "MLN", "PP")
treatment_vec <- c("UT", "EAE")
for (tissue in tissue_vec) {
  cells_to_use <- so_all$tissue == tissue
  # prefilter genes
  gene_prefilter <- 0.1 # a gene should be detected by this proportion in at least one of the tissues
  prop_mat <- sapply(treatment_vec, function(x){
    Matrix::rowMeans(so_all[["RNA"]]@counts[,cells_to_use & so_all$treatment == x] > 0) 
  }) 
  genes_to_use <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
  
  #### run edgeR
  meta_sub <- so_all@meta.data[cells_to_use,c("treatment", "batch", "nFeature_RNA")]
  # Create dummy variables
  meta_sub$treatment <- factor(meta_sub$treatment, levels = treatment_vec)
  meta_sub$batch <- factor(meta_sub$batch)
  contrasts(meta_sub$batch) = contr.sum(nlevels(meta_sub$batch))
  # Create a DGEList data object.
  dgeFull <- DGEList(so_all[['RNA']]@counts[genes_to_use, cells_to_use])
  # Estimate the normalization factors
  dgeFull <- calcNormFactors(dgeFull, method="TMM")
  # compute detection rate [here I use scaled number of detected genes]
  meta_sub$cdr <- scale(meta_sub$nFeature_RNA)
  # create design matrix
  design <- model.matrix(~ treatment + batch + cdr, data = meta_sub)
  # Estimate dispersion
  dgeFull <- estimateDisp(dgeFull, design = design)
  # Perform QLF tests
  fit <- glmQLFit(dgeFull, design = design)
  qlf <- glmQLFTest(fit, coef = 2)
  results <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
  saveRDS(results, file = paste0(dir_out, tissue, "_results.rds"))
}


