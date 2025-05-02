# Load required packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) BiocManager::install("TreeSummarizedExperiment")
if (!requireNamespace("mia", quietly = TRUE)) BiocManager::install("mia")

# Load libraries
library(ANCOMBC)
library(tidyverse)
library(TreeSummarizedExperiment)
library(mia)
library(readr)
library(ggrepel)

# Set working directory
setwd("~/Downloads/ANSC 516/Final Project_Apekshya/data/core-metrics-results")

# Create output directory if needed
if (!dir.exists("output/ANCOM_BC")) dir.create("output/ANCOM_BC", recursive = TRUE)

# Load genus and metadata RDS files
genus <- readRDS("genus_table.rds")
metadata <- readRDS("metadata_table.rds")

# Clean sample names
colnames(genus) <- make.names(colnames(genus), unique = TRUE)
rownames(metadata) <- make.names(rownames(metadata), unique = TRUE)

# Align and validate metadata
metadata <- metadata[colnames(genus), , drop = FALSE]
stopifnot(all(colnames(genus) == rownames(metadata)))

# Build TreeSummarizedExperiment with explicit names
counts <- as.matrix(genus)
sample_ids <- colnames(counts)
colnames(counts) <- sample_ids
rownames(metadata) <- sample_ids

# Ensure names are explicitly set in colData
colData_df <- S4Vectors::DataFrame(metadata)
rownames(colData_df) <- sample_ids

# Build the TSE
tse <- TreeSummarizedExperiment(
  assays = SimpleList(counts = counts),
  colData = colData_df
)

colData_df <- S4Vectors::DataFrame(metadata)
rownames(colData_df) <- colnames(genus)

which(is.na(rownames(colData_df)))
anyDuplicated(rownames(colData_df))

colnames(genus)[1:5]
rownames(colData_df)[1:5]


# Run ANCOM-BC v2
ancom_res <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level = "none",
  fix_formula = "Intervention_Status",
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05
)

# Extract results
res_df <- ancom_res$res$diff_abn
lfc_df <- ancom_res$res$lfc
padj_df <- ancom_res$res$p_val_adj

# Prepare plot data
plot_df <- res_df %>%
  rownames_to_column("taxon") %>%
  mutate(lfc = lfc_df$`Intervention_StatusMedA - Post`,
         padj = padj_df$`Intervention_StatusMedA - Post`) %>%
  filter(!is.na(lfc) & !is.na(padj))

# Save results
write.csv(plot_df, "output/ANCOM_BC/ancombc_results.csv", row.names = FALSE)

# Volcano plot with labels
ggplot(plot_df, aes(x = lfc, y = -log10(padj), color = diff_abn)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_n(plot_df, 5, wt = abs(lfc)),
                  aes(label = taxon), size = 3, max.overlaps = 10) +
  labs(title = "Volcano Plot: ANCOM-BC",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted p-value") +
  theme_minimal() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  ggsave("output/ANCOM_BC/volcano_plot_labeled.pdf", width = 7, height = 5)

# Top 20 barplot
top20 <- plot_df %>%
  filter(diff_abn == TRUE) %>%
  arrange(desc(abs(lfc))) %>%
  head(20)

ggplot(top20, aes(x = reorder(taxon, lfc), y = lfc, fill = lfc > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Differentially Abundant Taxa",
       x = "Taxa",
       y = "Log2 Fold Change") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "darkorange")) +
  theme(legend.position = "none") +
  ggsave("output/ANCOM_BC/top20_barplot.pdf", width = 7, height = 5)
