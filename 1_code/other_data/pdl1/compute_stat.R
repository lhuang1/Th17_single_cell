## compute statistics to show that Pdl1 (Cd274) is a good marker
rm(list = ls())

setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
# library(pROC)
library(RColorBrewer)
library(epitools)

today <- "2020-07-30"
dir_out <- paste0("3_output/FigureS4/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

## load data
data_names <- c("Tibbitt", "Hemmers", "Arazi")
df_cd274 <- sapply(data_names, simplify = FALSE, FUN = function(i) {
  readRDS(paste0("2_pipeline/Pdl1_", i, "/df_cd274_2020-07-23.rds"))
})

## density plot
df_cd274_combined <- lapply(names(df_cd274), function(i) {
  x <- df_cd274[[i]]
  x$data <- i
  x
}) %>% Reduce(rbind, .)


## Odds ratios
OR <- sapply(data_names, function(x) {
  
  ORtable <- table(df_cd274[[x]]$expr > 0, df_cd274[[x]]$in_cluster) 
  OR_res <- oddsratio.wald(ORtable)
  
  c(OR_res$measure[2,], OR_res$p.value[2,])
})

plt_df <- data.frame(t(OR)) %>% rownames_to_column(var = "Data") %>% 
  mutate(log_estimate = log(estimate),
         log_lower = log(lower),
         log_upper = log(upper)) %>% 
  mutate(Data = factor(Data, levels = data_names))
  
p_OR <- ggplot(plt_df) +
  geom_point(aes(x = Data, y = estimate), size = 2) +
  geom_errorbar(aes(x = Data, ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylim(0, max(plt_df$upper)) +
  labs(x = "", y = "Odds Ratio") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

p_logOR <- ggplot(plt_df) +
  geom_point(aes(x = Data, y = log_estimate), size = 2) +
  geom_errorbar(aes(x = Data, ymin = log_lower, ymax = log_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylim(0, max(plt_df$log_upper) + 1) +
  labs(x = "", y = "Log Odds Ratio") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none")

pdf(paste0(dir_out, "Pdl1_non_zero_OR_", today, ".pdf"), width = 4, height = 2)
p_OR
p_logOR
dev.off()

