### Project-wide utility functions ####
dir_proj <- "/singerlab/linglin/Th17_single_cell/"

## create mapping between human and mouse homology
create_hm_mapping <- function() {
  library(dplyr)
  f_out <- paste0(dir_proj, "output/results/other/hm_mapping.rds")
  if (!file.exists(f_out)) {
    hm_mapping <- read.table(paste0(dir_proj, "data/gene_lists/HOM_MouseHumanSequence_20200122.rpt.txt"), 
                             header = T, sep = "\t", stringsAsFactors = F)
    human <- hm_mapping %>% filter(NCBI.Taxon.ID == 9606)
    mouse <- hm_mapping %>% filter(NCBI.Taxon.ID == 10090)
    
    h2m <- mouse$Symbol[match(human$HomoloGene.ID, mouse$HomoloGene.ID)]
    names(h2m) <- human$Symbol
    
    m2h <- human$Symbol[match(mouse$HomoloGene.ID, human$HomoloGene.ID)]
    names(m2h) <- mouse$Symbol
    
    saveRDS(list("h2m" = h2m, "m2h" = m2h), file = f_out)
  }
  return(f_out)
}

## human to mouse gene names
human_to_mouse <- function(genes) {
  mh_mapping <- readRDS(create_hm_mapping())
  if (is.atomic(genes)) {  ## "genes" is one vector
    mapped <- mh_mapping$h2m[genes]
    mapped <- as.character(mapped[!is.na(mapped)])
  } else { ## "genes" is a list of multiple gene vectors
    mapped <- lapply(genes, function(x){
      tmp <- mh_mapping$h2m[x]
      return(as.character(tmp[!is.na(tmp)]))
    })
  }
  return(mapped)
}

## mouse to human gene names
mouse_to_human <- function(genes) {
  mh_mapping <- readRDS(create_hm_mapping())
  names(mh_mapping$m2h) <- toupper(names(mh_mapping$m2h))
  if (is.atomic(genes)) {  ## "genes" is one vector
    mapped <- mh_mapping$m2h[toupper(genes)]
    mapped <- as.character(mapped[!is.na(mapped)])
  } else { ## "genes" is a list of multiple gene vectors
    mapped <- lapply(genes, function(x){
      tmp <- mh_mapping$m2h[toupper(x)]
      return(as.character(tmp[!is.na(tmp)]))
    })
  }
  return(mapped)
}