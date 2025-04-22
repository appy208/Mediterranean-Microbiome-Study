
# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(dplyr)
library(ggpubr)
library(tidyr)


#Load qiime2R only after loading ggpubr
library(qiime2R)

setwd("~/Downloads/ANSC 516/Final Project_Apekshya/data")
getwd()

# Read the valid metadata file
library(readxl)
Final_Metadata_With_Intervention_Status <- read_excel("Final_Metadata_With_Intervention_Status.xlsx")
View(Final_Metadata_With_Intervention_Status)

#Check names of all columns in the metadata
colnames(metadata)

#Check for hidden characters
dput(colnames(metadata))

# Load cleaned metadata file
metadata <- Final_Metadata_With_Intervention_Status

#Renaming the column in metadata to SampleID to match codes
metadata <- metadata %>% rename(SampleID = `#SampleID`)
head(metadata)
str(metadata)

# Load alpha diversity metrics
shannon <- read_qza("core-metrics-results/shannon_vector.qza")$data
faith <- read_qza("core-metrics-results/faith_pd_vector.qza")$data
observed <- read_qza("core-metrics-results/observed_features_vector.qza")$data
evenness <- read_qza("core-metrics-results/evenness_vector.qza")$data

colnames(shannon)
colnames(faith)
colnames(observed)
colnames(evenness)

library(tibble)  # in case it's not loaded
shannon <- rownames_to_column(shannon, var = "SampleID")
faith <- rownames_to_column(faith, var = "SampleID")
observed <- rownames_to_column(observed, var = "SampleID")
evenness <- rownames_to_column(evenness, var = "SampleID")

colnames(shannon)
head(shannon$SampleID)
head(metadata$SampleID)
str(metadata)

colnames(faith)
head(faith$SampleID)

colnames(observed)
head(observed$SampleID)

colnames(evenness)
head(evenness$SampleID)

# Merge all alpha metrics into a single dataframe
alpha_df <- metadata %>%
  left_join(shannon, by = "SampleID") %>%
  left_join(faith, by = "SampleID") %>%
  left_join(observed, by = "SampleID") %>%
  left_join(evenness, by = "SampleID")

# Pivot to long format for ggplot
alpha_df_long <- alpha_df %>%
  pivot_longer(cols = c(shannon_entropy, faith_pd, observed_features, pielou_evenness),
               names_to = "Metric", values_to = "Value")
library(ggpubr)

# Define metrics and titles
alpha_metrics <- c("shannon_entropy", "faith_pd", "observed_features", "pielou_evenness")
titles <- c("Shannon Diversity", "Faith's Phylogenetic Diversity", "Observed Features", "Pielou's Evenness")


# Generate each plot with Kruskal Walis Statistics and collect in a list
alpha_plots <- list()

for (i in seq_along(alpha_metrics)) {
  p <- ggplot(alpha_df, aes(x = Intervention_Status, y = .data[[alpha_metrics[i]]], fill = Intervention_Status)) +
    geom_boxplot() +
    stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = "top") +
    theme_minimal() +
    labs(title = titles[i], y = titles[i], x = "") +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  alpha_plots[[i]] <- p
}

# Arrange in a 2x2 grid
combined_plot <- ggarrange(plotlist = alpha_plots, ncol = 2, nrow = 2, 
                           common.legend = TRUE, legend = "bottom")

# Save to PDF or PNG
ggsave("output/AlphaDiversity_Combined.pdf", combined_plot, width = 10, height = 8)
ggsave("output/AlphaDiversity_Combined.png", combined_plot, width = 10, height = 8, dpi = 300)


# Plot alpha diversity
#ggplot(alpha_df_long, aes(x = Intervention_Status, y = Value, fill = Intervention_Status)) +
#  geom_boxplot() +
 # facet_wrap(~ Metric, scales = "free_y") +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#  labs(title = "Alpha Diversity Metrics", y = "Diversity Value", x = "Intervention Group")


#alpha_metrics <- c("shannon_entropy", "faith_pd", "observed_features", "pielou_evenness")
#titles <- c("Shannon Diversity", "Faith's PD", "Observed Features", "Pielou's Evenness")

# Set output file and dimensions
#pdf("AlphaDiversity_Boxplots.pdf", width = 8, height = 6)

# Generate and save each plot
#for (i in seq_along(alpha_metrics)) {
 # p <- ggplot(alpha_df, aes(x = Intervention_Status, y = .data[[alpha_metrics[i]]], fill = Intervention_Status)) +
  #  geom_boxplot() +
   # stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = "top") +
    #theme_minimal() +
 #   theme_minimal() +
  #  labs(title = titles[i], y = titles[i], x = "") +
   # scale_fill_brewer(palette = "Set2") +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 # print(p)
#}

#Beta-Diversity Analysis & Plots
# Load necessary libraries
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tibble)

metadata <- Final_Metadata_With_Intervention_Status %>%
  rename(SampleID = `#SampleID`)

# File names and plot titles #Loads .qza PCoA files
ordination_files <- c("jaccard_pcoa_results.qza", 
                      "bray_curtis_pcoa_results.qza",
                      "unweighted_unifrac_pcoa_results.qza", 
                      "weighted_unifrac_pcoa_results.qza")

ordination_titles <- c("Jaccard", "Bray-Curtis", "Unweighted UniFrac", "Weighted UniFrac")

# Open a PDF device to save plots
pdf("BetaDiversity_PCoA_Plots.pdf", width = 8, height = 6)

# Loop through each ordination and plot
for (i in seq_along(ordination_files)) {
  
# Load QZA and extract data
  ord <- read_qza(paste0("core-metrics-results/", ordination_files[i]))$data
  
# Joins ordination files with cleaned metadata
  ord_df <- ord$Vectors %>%
    left_join(metadata, by = "SampleID")
  
  # Create plot
  p <- ggplot(ord_df, aes(x = PC1, y = PC2, color = Intervention_Status)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste0(ordination_titles[i], " PCoA"),
      x = paste0("PC1 (", round(ord$ProportionExplained[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(ord$ProportionExplained[2] * 100, 1), "%)")
    ) +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print plot to PDF
  print(p)
}

dev.off()

#PERMANOVA For Beta-diversity

library(qiime2R)
library(vegan)
library(dplyr)
library(tibble)

# Load and prepare metadata
metadata <- Final_Metadata_With_Intervention_Status %>%
  rename(SampleID = `#SampleID`) %>%
  column_to_rownames("SampleID")

# Define distance matrices and their names
distance_files <- c("bray_curtis_distance_matrix.qza",
                    "jaccard_distance_matrix.qza",
                    "unweighted_unifrac_distance_matrix.qza",
                    "weighted_unifrac_distance_matrix.qza")

distance_names <- c("Bray-Curtis", "Jaccard", "Unweighted UniFrac", "Weighted UniFrac")

# Run PERMANOVA for each
for (i in seq_along(distance_files)) {
  
  # Load distance matrix
  dist_matrix <- read_qza(paste0("core-metrics-results/", distance_files[i]))$data
  dist_mat <- as.matrix(dist_matrix)
  
  # Match samples
  common_ids <- intersect(rownames(dist_mat), rownames(metadata))
  dist_mat <- dist_mat[common_ids, common_ids]
  meta_filtered <- metadata[common_ids, , drop = FALSE]
  
  # Run PERMANOVA
  cat("\n--- PERMANOVA for", distance_names[i], "---\n")
  result <- adonis2(as.dist(dist_mat) ~ Intervention_Status, data = meta_filtered)
  print(result)
}

