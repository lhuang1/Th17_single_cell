setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(reticulate)

np <- import("numpy")
# data reading
mat <- np$load("0_data/other/Schafflick_GSE130119_RAW/batchid.npy", allow_pickle = T)
mat

library(Seurat)
so = readRDS("2_pipeline/human_shafflick/2020-02-26/so_combined.rds")

dim(so)
table(so$tissue)
table(so$orig.ident)


p_vec <- c('MS19270', 'MS49131','MS71658','MS60249','MS74594', 'PTC32190','PTC41540','PTC85037','PST83775','PST95809')

p_vec_my <- unique(so$individual[so$tissue == "PBMCs"])

setdiff(p_vec, p_vec_my)
setdiff(p_vec_my, p_vec)


library(openxlsx)
d = read.xlsx("0_data/other/Schafflick_GSE130119_RAW/41467_2019_14118_MOESM2_ESM.xlsx")
