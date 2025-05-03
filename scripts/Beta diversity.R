# Load required libraries
library(tidyverse)
library(vegan)
library(qiime2R)
library(ggplot2)

# Set working directory to your core-metrics-results folder
setwd("~/Downloads/ANSC 516/Final Project_Apekshya/data/core-metrics-results")

# Load metadata
metadata <- read_tsv("Final_Metadata_With_Intervention_Status.tsv")

# Set rownames
row.names(metadata) <- metadata$SampleID

# Create output directory
if (!dir.exists("beta_plots")) dir.create("beta_plots")

# Read PCoA results
bc_PCoA <- read_qza("bray_curtis_pcoa_results.qza")

# Merge with metadata
bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = "SampleID")

colnames(bc_meta)

# Define color palette
intervention_colors <- c("Baseline" = "black", "MedA - Post" = "blue", "MedWL - Post" = "green")

# Plot with ellipses
ggplot(bc_meta, aes(x = PC1, y = PC2, color = .data$Intervention_Status)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "t", level = 0.95) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(100 * bc_PCoA$data$ProportionExplained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(100 * bc_PCoA$data$ProportionExplained[2], 2), "%)")) +
  scale_color_manual(values = intervention_colors, name = "Intervention Status")

# Save plot
ggsave("beta_plots/BrayCurtis_InterventionStatus_PCoA.pdf", width = 7, height = 6)
ggsave("beta_plots/BrayCurtis_InterventionStatus_PCoA.png", width = 7, height = 6, dpi = 300)
