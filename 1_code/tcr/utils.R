meta_to_clonotype <- function(meta) {
  treatment_vec <- c("UT", "EAE")
  batch_vec <- c("b6", "b7", "b8")
  mouse_vec <- c("m1", "m2")
  meta_clonotype <- data.frame(matrix(NA, nrow = sum(!is.na(meta$raw_clonotype_id)), ncol = 2))
  rownames(meta_clonotype) <- rownames(meta)[!is.na(meta$raw_clonotype_id)]
  colnames(meta_clonotype) <- c("clonotype_id", "clone_size")
  for (treatment in treatment_vec) {
    for (batch in batch_vec) {
      for (mouse in mouse_vec) {
        meta_sub <- meta[meta$treatment == treatment & meta$batch == batch & meta$mouse == mouse & !is.na(meta$raw_clonotype_id),]
        if (nrow(meta_sub) == 0) next ## no such treatment, batch, mouse combination
        cat(treatment, batch, mouse, "\n")
        clonotypes <- meta_sub %>% rownames_to_column(var = "cellname") %>% 
          group_by(TRA, TRB) %>% summarize(cells_in_clone = paste(unique(cellname), collapse = ";"), clone_size = length(unique(cellname))) %>% 
          arrange(-clone_size) %>%
          data.frame()
        clonotypes$clonotype_id <- paste("clonotype", seq_along(clonotypes[,1]), sep = "")
        cells_list <- strsplit(clonotypes$cells_in_clone, split = ";")
        # tmp <- lapply(seq_along(clonotypes$clonotype_id), function(i) {
        #   data.frame(cell = cells_list[[i]], clonotype = clonotypes$clonotype[i])
        # }) %>% Reduce(rbind, .) %>% mutate(treatment = treatment, batch = batch, mouse = mouse)
        # meta_clonotype <- rbind(meta_clonotype, tmp)
        for (i in seq_along(clonotypes$clonotype_id)) {
          meta_clonotype[cells_list[[i]],"clonotype_id"] <-  clonotypes$clonotype_id[i]
          meta_clonotype[cells_list[[i]],"clone_size"] <-  clonotypes$clone_size[i]
        }
      }
    }
  }
  return(meta_clonotype)
}




## compute and compare pct cells sharing in clusters
compare_pcts <- function(tissue_x, cluster_x, tissue_y, cluster_y, ann_type) {
  meta_x <- meta_noNA %>% filter(tissue == tissue_x) %>% 
    mutate(in_cluster = tissue_cluster %in% cluster_x)
  meta_y <- meta_noNA %>% filter(tissue == tissue_y) %>% 
    mutate(in_cluster = tissue_cluster %in% cluster_y)
  
  sharing_info <- list()
  for (s in levels(meta_noNA$sample)) {
    meta_x_sub <- meta_x %>% filter(sample == s)
    meta_y_sub <- meta_y %>% filter(sample == s)
    if (nrow(meta_y_sub) == 0 | nrow(meta_x_sub) == 0) next()
    clones_in_cluster_y <- meta_y_sub$clonotype_id[meta_y_sub$in_cluster] %>% as.character() %>% unique()
    clones_out_cluster_y <- meta_y_sub$clonotype_id[!meta_y_sub$in_cluster] %>% as.character() %>% unique()
    x_in_y <- table(shared = meta_x_sub$clonotype_id %in% clones_in_cluster_y, x_in_cluster = meta_x_sub$in_cluster) %>% as.data.frame() %>% 
      group_by(x_in_cluster) %>% mutate(Total = sum(Freq)) %>% ungroup() %>% mutate(Frac = Freq / Total)
    x_in_y$y_in_cluster <- TRUE
    x_out_y <- table(shared = meta_x_sub$clonotype_id %in% clones_out_cluster_y, x_in_cluster = meta_x_sub$in_cluster) %>% as.data.frame() %>% 
      group_by(x_in_cluster) %>% mutate(Total = sum(Freq)) %>% ungroup() %>% mutate(Frac = Freq / Total)
    x_out_y$y_in_cluster <- FALSE
    sharing_info[[s]] <- cbind(rbind(x_in_y, x_out_y), sample = s)
  }
  ## pooled
  clones_in_cluster_y <- meta_y$pooled_clonotype_id[meta_y$in_cluster] %>% as.character() %>% unique()
  clones_out_cluster_y <- meta_y$pooled_clonotype_id[!meta_y$in_cluster] %>% as.character() %>% unique()
  x_in_y <- table(shared = meta_x$pooled_clonotype_id %in% clones_in_cluster_y, x_in_cluster = meta_x$in_cluster) %>% as.data.frame() %>% 
    group_by(x_in_cluster) %>% mutate(Total = sum(Freq)) %>% ungroup() %>% mutate(Frac = Freq / Total)
  x_in_y$y_in_cluster <- TRUE
  x_out_y <- table(shared = meta_x$pooled_clonotype_id %in% clones_out_cluster_y, x_in_cluster = meta_x$in_cluster) %>% as.data.frame() %>% 
    group_by(x_in_cluster) %>% mutate(Total = sum(Freq)) %>% ungroup() %>% mutate(Frac = Freq / Total)
  x_out_y$y_in_cluster <- FALSE
  sharing_info[["pooled"]] <- cbind(rbind(x_in_y, x_out_y), sample = "pooled")
  
  
  sharing_info <- sharing_info %>% Reduce(rbind,.)
  
  plt_df <- sharing_info %>% filter(shared == "TRUE") %>% 
    mutate(x_in_cluster = paste0(tissue_x, "_", ifelse(x_in_cluster == "TRUE", ann_type, paste0("non-", ann_type))),
           y_in_cluster = paste0(tissue_y, "_", ifelse(y_in_cluster == "TRUE", ann_type, paste0("non-", ann_type))))
  
  p <- ggplot(plt_df, aes(x = x_in_cluster, y = Frac, group = sample, color = sample)) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~y_in_cluster) +
    labs(x = "", y = paste0("Fraction Shared with ", tissue_y), color = "Mouse") +
    ggtitle(paste0(ann_type, ": Fraction of ", tissue_x, " Shared with ", tissue_y)) +
    # geom_signif(comparisons = list(unique(plt_df$x_in_cluster)), map_signif_level = F, vjust = 0, color = "black", margin_top = 0.1) +
    theme_bw() +
    theme(panel.border = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.line = element_line(),
          axis.text.x = element_text(size = 10))
  print(p)
  return(sharing_info)
}
