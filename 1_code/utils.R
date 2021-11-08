#### utility functions ####

## create mapping between human and mouse homology
create_hm_mapping <- function() {
  library(dplyr)
  f_out <- "/singerlab/linglin/Th17_single_cell_eae_ut/2_pipeline/other/hm_mapping.rds"
  if (!file.exists(f_out)) {
    hm_mapping <- read.table("/singerlab/linglin/Th17_single_cell_eae_ut/0_data/gene_lists/HOM_MouseHumanSequence_20200122.rpt.txt", 
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


## tissue names mapping
tissue_informal_names <- function(tissues) {
  informal_names <- c("SPL" = "spleen", 
                      "PP" = "PP", 
                      "MLN" = "mLN", 
                      "SI" = "SI", 
                      "COL" = "colon", 
                      "CNS" = "CNS",
                      "dLN" = "DLN")
  return(informal_names[tissues])
}


## save and plot marker genes 
save_and_plot_markers <- function(so, so_markers, dir_markers, pt.size = 0.2){
  cluster_names <- unique(so_markers$cluster)
  for (x in cluster_names) {
    cat("Cluster", x, "...\n")
    cluster_markers = so_markers %>% filter(cluster == x) #%>% arrange(desc(avg_logFC))
    write.table(cluster_markers, 
                file = paste0(dir_markers, "/FILES/Cluster_Markers_", x, ".txt"), quote = F)
    cluster_markers_no_rps_rpl  = cluster_markers %>% 
      filter(!(grepl("^Rps", gene) | grepl("^Rpl", gene)))
    write.table(cluster_markers_no_rps_rpl, 
                file = paste0(dir_markers, "/FILES/Cluster_Markers_", x, "_no_rps_rpl.txt"), quote = F)
    cluster_markers_thresh03  = cluster_markers %>% 
      filter(pct.2 < 0.3)
    write.table(cluster_markers_thresh03, 
                file = paste0(dir_markers, "/FILES/Cluster_Markers_", x, "_thresh03.txt"), quote = F)
    to_plot <- cluster_markers$gene[1:min(18, nrow(cluster_markers))]
    pdf(file = paste0(dir_markers, "/PLOTS/Cluster_Markers_", so@active.assay, "_", x, ".pdf"))
    for (to_plot_i in to_plot) {
      print(VlnPlot(object = so, features = to_plot_i, pt.size = pt.size))
      print(FeaturePlot(object = so, features = to_plot_i, 
                        min.cutoff = "q5", 
                        max.cutoff = "q95",
                        reduction = "umap", pt.size = pt.size))
    }
    dev.off()
  }
  return(NULL)
}


SPL_EAE_cluster_0_1_colors <- function(x){
  colors <- c("SPL_0" = "#1E90FF", "SPL_1" = "#DC143C")
  return(colors[x])
  
}

SPL_EAE_cluster_colors <- function(x){
  colors <- c("Cxcr6" = "#DC143C", "Slamf6" = "#1E90FF")
  return(colors[x])

}

SPL_EAE_all_cluster_colors <- function(x = NULL){
  colors <- c("0" = "#1E90FF", "1" = "#DC143C", "2" = "#FFA500", "3" = "#9ACD32", "4" = "#BA55D3")
  if (is.null(x)) return(colors)
  return(colors[x])
}


tissue_colors <- function(tissues = NULL){
  color_vec <- c(
    "SPL" = "#E76BF3", # pink FF99CC(lighter pink)
    "MLN" = "#A3A500", #pickle
    "PP" = "#00BF7D", # green
    "SI" = "#00B0F6", # blue
    "COL" = "#F8766D", # red f65348
    "CNS" = "#FFA500", # orange 
    "DLN" = "#6495ED"
  )
  if (is.null(tissues)) {
    return(color_vec)
  } else {
    return(color_vec[tissues])
  }
}

gfp_colors <- function() {
  color_vec <- c(
    "Ex Th17" = "gray70",
    "Current Th17" = "#006600"
  )
  return(color_vec)
}

ggplot_default_color_descrete <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col_corr <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))
col_corr_no_white <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582" , "#FDDBC7", "#D1E5F0", "#92C5DE",
                                        "#4393C3", "#2166AC", "#053061"))
col_corr_pos <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF"))
col_corr_pos_no_white <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7"))
col_corr_pos_skip <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))
col_corr_neg <- colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))

