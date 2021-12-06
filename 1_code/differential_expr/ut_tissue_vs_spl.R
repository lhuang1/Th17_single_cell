################################################################
######## Differential Expression: tissue vs. SPL ############
################################################################

## Identify differentially expressed genes comparing each tissue to SPL (UT, GFPall)
## edgeR/QLF with cellular detection rate as covariate
## perform test for each tissue-batch vs. SPL all-batch-average

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

library(Seurat)
library(dplyr)
library(edgeR)
set.seed(1)

## configuration
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()
  trt_gfp <- "UT_GFPpos" # or "UT_GFPall", "EAE_UT_GFPall"; treatment and GFP status combination
  prep_date <- "2020-03-24"
} else {
  today <- cargs[1]
  trt_gfp <- cargs[2]
  prep_date <- cargs[3]
}
dir_out <- paste0("2_pipeline/differential_expression/", trt_gfp, "_vs_SPL/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_dominant_TCR_", prep_date, ".rds"))

if (trt_gfp == "UT_GFPpos") {
  keep_cells_sub <- colnames(so_all)[so_all$treatment == "UT" & so_all$batch %in% paste0("b", 1:5) & so_all$GFP_positive]
} else if (trt_gfp == "UT_GFPall") {
  keep_cells_sub <- colnames(so_all)[so_all$treatment == "UT" & so_all$batch %in% paste0("b", 1:5)]
} else if (trt_gfp == "EAE_UT_GFPall") {
  keep_cells_sub <- colnames(so_all)
  so_all$tissue_treatment <- paste(so_all$tissue, so_all$treatment, sep = "_")
}
so <- so_all[,keep_cells_sub]
so@meta.data <- droplevels(so@meta.data)
rm(so_all)
DefaultAssay(object = so) <- "RNA"

########Differential Expression#########
## pre-filter genes
gene_prefilter <- 0.1 # a gene should be detected by this proportion in at least one of the tissues
if (trt_gfp == "UT_GFPpos") {
  prop_mat <- sapply(unique(so$tissue), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue == x] > 0) 
  }) 
  keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
  so <- subset(so, features = which(keep))
} else if (trt_gfp == "UT_GFPall") {
  prop_mat <- sapply(unique(so$tissue), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue == x] > 0) 
  }) 
  keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
  ## make sure genes used in GFP pos are also included so can compare results
  prop_mat_gfppos <- sapply(unique(so$tissue), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue == x & so$GFP_positive] > 0) 
  }) 
  keep_gfppos <- (prop_mat_gfppos > gene_prefilter) %>% Matrix::rowSums(.) > 0
  so <- subset(so, features = which(keep | keep_gfppos))
} else if (trt_gfp == "EAE_UT_GFPall") {
  prop_mat <- sapply(unique(so$tissue_treatment), function(x){
    Matrix::rowMeans(so@assays$RNA@counts[,so$tissue_treatment == x] > 0) 
  }) 
  keep <- (prop_mat > gene_prefilter) %>% Matrix::rowSums(.) > 0
  so <- subset(so, features = which(keep))
}

## run edgeR
so@meta.data$batch <- factor(so@meta.data$batch, levels = paste0("b", 1:5)) %>% droplevels()
so@meta.data$tissue <- factor(so@meta.data$tissue, levels = c("SPL", "MLN", "PP", "SI", "COL")) %>% droplevels()
so@meta.data$tissue_batch <- factor(paste(so@meta.data$tissue, so@meta.data$batch, sep = "_"))
# Create dummy variables
contrasts(so@meta.data$batch) = contr.sum(nlevels(so@meta.data$batch))
# Create a DGEList data object.
dgeFull <- DGEList(so@assays$RNA@counts, group=so$tissue, samples = so$orig.ident)
# Estimate the normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")
# compute detection rate [here I use scaled number of detected genes from the original unfiltered data]
so@meta.data$cdr <- scale(so$nFeature_RNA)
# create design matrix
design <- model.matrix(~ tissue + batch + cdr, data = so@meta.data)
# Estimate dispersion
dgeFull <- estimateDisp(dgeFull, design = design)
# Perform QLF tests
fit <- glmQLFit(dgeFull, design = design)
# Construct contrast matrix
contrast_mat <- matrix(0, nrow = (nlevels(so$tissue_batch) + nlevels(so$tissue) - 1), ncol = ncol(design),
                       dimnames = list(c(levels(so$tissue_batch), setdiff(levels(so$tissue), "SPL")), 
                                       colnames(design)))
for (tb in rownames(contrast_mat)) {
  tb_vec <- unlist(strsplit(tb, split = "_"))
  if (length(tb_vec) == 2) { # tissue and batch
    if (tb_vec[1] != "SPL") {
      if (tb_vec[2] != "b5") {
        coef_t <- paste0("tissue", tb_vec[1])
        coef_b <- paste0("batch", sub("b", "", tb_vec[2]))
        contrast_mat[tb, c(coef_t, coef_b)] <- 1
      } else {
        coef_t <- paste0("tissue", tb_vec[1])
        contrast_mat[tb, coef_t] <- 1
        coef_b <- paste0("batch", 1:4)
        contrast_mat[tb, coef_b] <- -1
      }
    } else { # other tissues
      if (tb_vec[2] != "b5") {
        coef_b <- paste0("batch", sub("b", "", tb_vec[2]))
        contrast_mat[tb, coef_b] <- 1
      } else {
        coef_b <- paste0("batch", 1:4)
        contrast_mat[tb, coef_b] <- -1
      }
    }
  } else if (length(tb_vec) == 1) { # tissue
    coef_t <- paste0("tissue", tb_vec[1])
    contrast_mat[tb, coef_t] <- 1
  }
}

results <- lapply(rownames(contrast_mat), function(tb){
  cat(tb, "...\n")
  qlf <- glmQLFTest(fit, contrast = contrast_mat[tb,])
  # Extract summaries of differential expression statistics.
  res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
  return(res$table)
})
names(results) <- rownames(contrast_mat)
saveRDS(results, file = paste0(dir_out, "results.rds"))