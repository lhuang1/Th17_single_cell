###########################################
######## clean TCR information ############

## Chain and cell QC
# 1). Filter chains by "cell" quality, productivity, confidence and chain type
# 2). Identify dominant chains (for alpha and beta respectively) based on UMI counts: 
#   (1). only one chain: mark this chain as dominant
# (2). only two chains: mark the more abundant chain as dominant, or both as dominant if equal UMIs.
# (3). more than two chain, and the UMI count of the most abundant one is >=2 times the second abundant: mark the most abundant chain as dominant (if >=2 chains are the most abundant, no dominant chain.)
# (4). more than two chains, and the UMI count of the most abundant one is <2 time the second abundant: no dominant chain
# 3). Filter cells by number of dominant chains: 
#   (1) remove cells with no dominant beta chains
# (2) remove cells with no dominant alpha chains (unless there were no chains meet QC criteria to begin with)

rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(dplyr)
set.seed(1)

#### configuration ####
today <- Sys.Date()
dir_out <- paste0("2_pipeline/preprocessing/tcr/", today, "/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

## load data
batch_vec <- paste0("batch", 6:9)
treatment_vec <- c("EAE", "UT")
mouse_vec <- c("m1", "m2")
for (batch in batch_vec) {
  for (treatment in treatment_vec) {
    for (mouse in mouse_vec) {
      sample_name <- paste(batch, "hashed", treatment, mouse, sep = "_")
      dir_vdj <- paste0(dir_proj, "data/single_cell/VDJ/", sample_name, "_all_contig_annotations.csv")
      if (!file.exists(dir_vdj)) next()
      cat("Processing", sample_name, "...\n")
      
      contigs_raw <- read.csv(dir_vdj, stringsAsFactors = F)
      
      ## transform cell names to match the format in expression data
      tmp <- strsplit(as.character(contigs_raw$barcode), split = '-')
      stopifnot(length(unique(sapply(tmp, function(x)x[2]))) == 1) # make sure all have same suffix, so it's safe to remove
      contigs_raw$barcode_nosuffix <- paste0(sapply(tmp, function(x) x[1]))
      
      ## filter out: is_cell == False; high_confidence == False; productive != True; chain not "TRA" or "TRB"
      contigs <- contigs_raw %>% filter(is_cell == "True" & high_confidence == "True" & productive == "True" & chain %in% c("TRA", "TRB")) %>% droplevels()
      
      ## identify dominant chains
      contigs <- contigs %>% group_by(barcode, chain) %>% 
        mutate(n_chains = n(), 
               sec_max_umi = sort(umis, decreasing = T)[2],
               is_dominant = ifelse(n_chains == 1, T, ifelse(n_chains == 2, umis == max(umis), (umis / sec_max_umi) >= 2)),
               n_dominant = sum(is_dominant))
      ## filter chains and cells
      contigs_dominant <- contigs %>% ungroup() %>% group_by(barcode) %>% 
        mutate(keep_cell = any(chain[is_dominant] == "TRB") & (any(chain[is_dominant] == "TRA") | all(chain != "TRA"))) %>% 
        ungroup() %>% filter(keep_cell, is_dominant) %>% droplevels()
      saveRDS(contigs_dominant, file = paste0(dir_out, sample_name, "_filtered_contigs_dominant.rds"))
    }
  }
}