#########################################################################
######## Differential Expression: spleen EAE cluster 0 vs. 1 ############
#########################################################################

## Identify differentially expressed genes comparing cluster 0 and 1 in SPL (EAE, GFP all)
## edgeR/QLF with cellular detection rate as covariate

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

source("1_code/utils.R")
source("1_code/differential_expr/utils.R")

library(Seurat)
library(dplyr)
library(tibble)
library(edgeR)
set.seed(1)

# configuration
today <- Sys.Date()
prep_date <- "2020-03-24"
clustering_date <- "2020-03-25"
dir_out <- paste0("2_pipeline/differential_expression/SPL_EAE_0_1/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_dominant_TCR_", prep_date, ".rds"))
meta <- read.table(paste0("2_pipeline/clustering/EAE_GFPall_SPL/", clustering_date, "/FILES/meta_data.csv"), sep = " ")
so <- so_all[,rownames(meta)[meta$seurat_clusters %in% c(0, 1)]]
so@meta.data <- droplevels(meta[colnames(so),])
rm(so_all)
DefaultAssay(object = so) <- "RNA"

########Differential Expression#########
so@meta.data$batch <- factor(so@meta.data$batch)
so@meta.data$seurat_clusters <- factor(so@meta.data$seurat_clusters)
# Create dummy variables
contrasts(so@meta.data$batch) = contr.sum(nlevels(so@meta.data$batch))
# Create a DGEList data object.
dgeFull <- DGEList(so[['RNA']]@counts)
# Estimate the normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")
# compute detection rate [here I use scaled number of detected genes]
so@meta.data$cdr <- scale(so$nFeature_RNA)
# create design matrix
design <- model.matrix(~ seurat_clusters + batch + cdr, data = so@meta.data)
# Estimate dispersion
dgeFull <- estimateDisp(dgeFull, design = design)
# Perform QLF tests
fit <- glmQLFit(dgeFull, design = design)
qlf <- glmQLFTest(fit, coef = 2)
results <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
saveRDS(results, file = paste0(dir_out, "results.rds"))

####### process results into signature ######
cxcr6 <- results$table %>% rownames_to_column(var = "gene") %>% filter(logFC < -log2(1.5) & FDR < 0.05) %>% arrange(FDR)
slamf6 <- results$table %>% rownames_to_column(var = "gene") %>% filter(logFC > log2(1.5) & FDR < 0.05) %>% arrange(FDR)
sig <- list(
  'Cxcr6_sc-plus' = cxcr6$gene,
  'Slamf6_sc-plus' = slamf6$gene
)
saveRDS(sig, file = paste0(dir_out, "sc_signature.rds"))
