library(Hmisc)
library(plyr)
library(reshape2)
library(qiime2R)
library(tidyverse)

setwd("~/Downloads/ANSC 516/Final Project_Apekshya/data")

if (!dir.exists("output/cooccurrence")) dir.create("output/cooccurrence", recursive = TRUE)

ASVs <- read_qza("core-metrics-results/rarefied_table.qza")
ASV_table <- as.data.frame(ASVs$data)

# Add ASV labels
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)]
ASV_table <- ASV_table[, -((ncol(ASV_table)-1):ncol(ASV_table))]

dataset <- as.data.frame(t(ASV_table))

# Load and join metadata
meta_raw <- read.delim("Final_Metadata_With_Intervention_Status.txt", 
                       header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)

# Assign correct column names manually
colnames(meta_raw) <- c("SampleID", "Sample_Name", "BMI", "BMI_Category", "Time_Point", 
                        "Randomization_Group", "Percent_Weight_Change", 
                        "MedDiet_Adherence_Change", "MediScore_Change",
                        "Classes_Attended", "Class_Attendance_Pct", "Age", 
                        "Intervention_Status")
metadata <- meta_raw  # use this going forward

# Create grouping variable
metadata$Intervention_ord <- factor(metadata$Intervention_Status, 
                                    levels = c("Baseline", "MedA - Post", "MedWL - Post"))

# Merge ASV dataset with metadata
dataset <- merge(metadata, dataset, by.x = "SampleID", by.y = 0)
datasetn <- dataset
datasetn[datasetn == 0] <- NA

my_column <- "Intervention_ord"
treatments <- unique(dataset[[my_column]])
num_metadata_columns <- ncol(metadata)

# Sample number thresholds per group (for correlation robustness)
n1 <- 8
n2 <- 8
n3 <- 8
q_cutoff <- 0.05
final_results <- data.frame()

for(i in 1:length(treatments)) {
  cat("Analyzing:", treatments[i], "\n")
  temp <- subset(dataset, get(my_column) == treatments[i])
  tempn <- subset(datasetn, get(my_column) == treatments[i])
  
  results <- rcorr(as.matrix(temp[, -(1:num_metadata_columns)]), type = "spearman")
  resultsn <- rcorr(as.matrix(tempn[, -(1:num_metadata_columns)]), type = "spearman")
  
  rhos <- results$r
  ps <- results$P
  ns <- resultsn$n
  
  ps_melt <- na.omit(melt(ps))
  ps_melt$qval <- p.adjust(ps_melt$value, method = "BH")
  colnames(ps_melt)[3] <- "pval"
  ps_sub <- subset(ps_melt, qval < q_cutoff)
  
  rhos_melt <- na.omit(melt(rhos)); colnames(rhos_melt)[3] <- "rho"
  ns_melt <- melt(ns); colnames(ns_melt)[3] <- "n"
  
  merged <- merge(ps_sub, rhos_melt, by = c("Var1", "Var2"))
  merged <- merge(merged, subset(ns_melt, n > get(paste0("n", i))), by = c("Var1", "Var2"))
  
  if(nrow(merged) > 0){
    merged$trt <- treatments[i]
    final_results <- rbind(final_results, merged)
  } else {
    cat("No significant correlations for:", treatments[i], "\n")
  }
}

# Load taxonomy
taxonomy <- read_qza("core-metrics-results/taxonomy.qza")
tax.clean <- parse_taxonomy(taxonomy$data)
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)) {
  for (j in 2:7) {
    if (tax.clean[i, j] == "") {
      tax.clean[i, j:7] <- paste0("unclassified_", tax.clean[i, j - 1])
      break
    }
  }
}

# Merge taxonomy
strong_results <- final_results
strong_results_taxa <- merge(strong_results, ASVkey, by.x = "Var1", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, ASVkey, by.x = "Var2", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.x", by.y = 0)
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.y", by.y = 0)

# Save all and per-group results
write.csv(strong_results_taxa, "output/cooccurrence/strong_results_taxa_all.csv", row.names = FALSE)
for (trt in treatments) {
  write.csv(
    subset(strong_results_taxa, trt == !!trt),
    file = paste0("output/cooccurrence/strong_results_taxa_", gsub(" ", "_", trt), ".csv"),
    row.names = FALSE
  )
}

# Example plot (replace ASV IDs as needed)
ggplot(subset(dataset, Intervention_ord == "Baseline"), aes(x = ASV48, y = ASV84)) + geom_point()

# Optional: create Cytoscape-ready edge table
cytoscape_edges <- strong_results_taxa %>%
  filter(abs(rho) >= 0.6 & qval < 0.05) %>%
  filter(Genus.x != Genus.y) %>%
  select(Source = Genus.x, Target = Genus.y, rho, qval, trt)

write.csv(cytoscape_edges, "output/cooccurrence/cytoscape_edges.csv", row.names = FALSE)

# Optional: create Cytoscape-ready node table
nodes <- unique(c(cytoscape_edges$Source, cytoscape_edges$Target))
cytoscape_nodes <- data.frame(ID = nodes)
write.csv(cytoscape_nodes, "output/cooccurrence/cytoscape_nodes.csv", row.names = FALSE)


# Compute Spearman correlation
corr <- cor(df$ASV15, df$ASV5, method = "spearman", use = "complete.obs")

# Create a labeled plot
library(ggplot2)

ggplot(df, aes(x = ASV48, y = ASV84)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed") +
  annotate("text", x = max(df$ASV48, na.rm = TRUE) * 0.8,
           y = max(df$ASV84, na.rm = TRUE) * 0.9,
           label = paste0("Spearman r = ", round(corr, 2)),
           size = 5, color = "black") +
  theme_minimal() +
  labs(
    title = "ASV48 vs ASV84",
    x = "ASV48 Abundance",
    y = "ASV84 Abundance"
  )


colnames(strong_results_taxa)

# Pick one pair to plot — change the row index to see different pairs
plot_row <- 1

# Extract genus names and ASV strings for selected pair
genus_x <- strong_results_taxa$Genus.x[plot_row]
genus_y <- strong_results_taxa$Genus.y[plot_row]
asv_x <- strong_results_taxa$ASVstring.x[plot_row]
asv_y <- strong_results_taxa$ASVstring.y[plot_row]

# Plot for a specific treatment group
ggplot(subset(dataset, Intervention_ord == "Baseline"), aes_string(x = asv_x, y = asv_y)) +
  geom_point() +
  labs(
    x = paste0("Genus: ", genus_x, " (", asv_x, ")"),
    y = paste0("Genus: ", genus_y, " (", asv_y, ")"),
    title = paste("Co-occurrence in Baseline Group")
  ) +
  theme_minimal()

library(ggplot2)
library(rlang)

####ASV84 VS ASV48
# STEP 2: Load dataset
cooc_data <- read_csv("output/cooccurrence/strong_results_taxa_all.csv")  # ✅ adjust path if needed

# STEP 3: Find rho for the exact ASV pair
target_pair <- subset(cooc_data, (Var1 == "ASV89" & Var2 == "ASV157") | (Var1 == "ASV157" & Var2 == "ASV89"))
print(target_pair$rho)
print(target_pair$trt)  # treatment group (important!)

# STEP 4: Get your merged abundance + metadata dataset
# Replace this with your own full merged dataframe (like 'dataset')
# It must contain ASV89, ASV157, and 'Intervention_ord' column
df <- dataset  # assuming you've already created this above

# STEP 5: Subset by treatment group from the `trt` column
df_subset <- subset(df, Intervention_ord == target_pair$trt[1])  # e.g., "MedWL - Post"

# STEP 6: Plot the exact same data as in correlation
asv_x <- "ASV89"
asv_y <- "ASV157"
genus_x <- target_pair$Genus.x[1]
genus_y <- target_pair$Genus.y[1]
rho_exact <- round(target_pair$rho[1], 2)

p <- ggplot(df_subset, aes(x = .data[[asv_x]], y = .data[[asv_y]])) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed") +
  annotate("text", x = max(df_subset[[asv_x]], na.rm = TRUE) * 0.8,
           y = max(df_subset[[asv_y]], na.rm = TRUE) * 0.9,
           label = paste0("Spearman r = ", rho_exact),
           size = 5, color = "black") +
  labs(
    x = paste0("Genus: ", genus_x, " (", asv_x, ")"),
    y = paste0("Genus: ", genus_y, " (", asv_y, ")"),
    title = paste0("Co-occurrence in ", target_pair$trt[1], " group")
  ) +
  theme_minimal()

print(p)
