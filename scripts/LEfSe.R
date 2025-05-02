# LEfSe-style Kruskal-Wallis Analysis with Effect Size on ASV-level data

# Load required libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)

# Set working directory
setwd("~/Downloads/ANSC 516/Final Project_Apekshya/data/core-metrics-results")

# Load ASV-level data, metadata, and taxonomy
asv <- readRDS("asv_table.rds")
metadata <- readRDS("metadata_table.rds")
taxonomy_table <- read_tsv("taxonomy.tsv")

# Clean sample names
colnames(asv) <- make.names(colnames(asv), unique = TRUE)
rownames(metadata) <- make.names(rownames(metadata), unique = TRUE)

# Ensure all sample IDs match
stopifnot(all(colnames(asv) %in% rownames(metadata)))

# Align metadata with ASV table
metadata <- metadata[colnames(asv), , drop = FALSE]

# Initialize LEfSe-style results
lefse_results <- lapply(rownames(asv), function(asv_id) {
  asv_id <- as.character(asv_id)
  values <- as.numeric(asv[asv_id, ])
  group <- metadata$Intervention_Status
  
  kw_result <- tryCatch(kruskal.test(values ~ group), error = function(e) NULL)
  if (!is.null(kw_result)) {
    kw_stat <- kw_result$statistic
    kw_df <- kw_result$parameter
    kw_p <- kw_result$p.value
  } else {
    kw_stat <- NA
    kw_df <- NA
    kw_p <- NA
  }
  
  eta2 <- function(x, g) {
    df <- data.frame(x = x, g = g)
    ss_total <- sum((x - mean(x))^2)
    ss_between <- sum(tapply(x, g, function(v) length(v) * (mean(v) - mean(x))^2))
    return(ss_between / ss_total)
  }
  eta <- tryCatch(eta2(values, group), error = function(e) NA)
  
  data.frame(
    asv_id = asv_id,
    chi_sq = kw_stat,
    df = kw_df,
    p_value = kw_p,
    effect_size = eta
  )
})

# Combine and adjust
lefse_df <- bind_rows(lefse_results) %>%
  mutate(q_value = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)

# Join with taxonomy
taxonomy_clean <- taxonomy_table %>% rename(asv_id = `Feature.ID`)
lefse_annotated <- lefse_df %>% left_join(taxonomy_clean, by = "asv_id")

# Save annotated results
if (!dir.exists("output")) dir.create("output")
write.csv(lefse_annotated, "output/lefse_asv_kruskal_annotated_results.csv", row.names = FALSE)

# Volcano plot with taxonomy label and color by significance
lefse_annotated <- lefse_annotated %>%
  mutate(significant = ifelse(q_value < 0.05, "Significant", "Not Significant"))

ggplot(lefse_annotated, aes(x = effect_size, y = -log10(q_value), color = significant)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = filter(lefse_annotated, q_value < 0.05 & effect_size > 0.01),
                  aes(label = Taxon), size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(title =         "LEfSe:Volcano Plot",
       x = "Effect Size (eta²)", y = "-log10(q-value)", color = "Significance") +
  theme_minimal()

ggsave("output/lefse_volcano_plot_asv_labeled.pdf", width = 7, height = 5)

# Condensed Kruskal-Wallis summary table for significant hits
lefse_summary <- lefse_annotated %>%
  filter(q_value < 0.05, effect_size > 0.01) %>%
  arrange(p_value) %>%
  mutate(rank = row_number()) %>%
  select(rank, asv_id, Taxon, chi_sq, df, p_value, q_value, effect_size)

write.csv(lefse_summary, "output/lefse_condensed_summary.csv", row.names = FALSE)

