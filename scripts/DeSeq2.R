# Updated full script with corrected metadata handling for taxa bar plots and DESeq2

# Load required libraries
library(phyloseq)
library(DESeq2)
library(qiime2R)
library(tidyverse)
library(readxl)

# Set working directory
setwd("~/Downloads/ANSC 516/Final Project_Apekshya/data")

# Create output folder
if(!dir.exists("output/taxa")) dir.create("output/taxa", recursive = TRUE)

# Read metadata from Excel and save as TSV
meta_df <- read_excel("Final_Metadata_With_Intervention_Status.xlsx")
write.table(meta_df, "Final_Metadata_With_Intervention_Status.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Import files into phyloseq
physeq <- qza_to_phyloseq(
  features = "core-metrics-results/rarefied_table.qza",
  tree = "rooted-tree.qza",
  taxonomy = "core-metrics-results/taxonomy.qza",
  metadata = "Final_Metadata_With_Intervention_Status.tsv"
)

# Extract tables
asv_table <- data.frame(otu_table(physeq), check.names = FALSE)
metadata <- data.frame(sample_data(physeq), check.names = FALSE)
taxonomy <- data.frame(tax_table(physeq), check.names = FALSE)

# Set factor levels
metadata$Intervention_ord <- factor(metadata$`MedA...Post`, levels = c("Baseline", "MedA - Post", "MedWL - Post"))
metadata$MedDt_ord <- factor(metadata$`X2..MedDt.A`, levels = c("1. MedDt-WL", "2. MedDt-A"))
metadata$Obese_ord <- factor(metadata$Obese, levels = c("Overweight", "Obese"))

# Push updated metadata back into phyloseq
sample_data(physeq)$Intervention_ord <- metadata$Intervention_ord
sample_data(physeq)$MedDt_ord <- metadata$MedDt_ord
sample_data(physeq)$Obese_ord <- metadata$Obese_ord

# Clean taxonomy table
tax.clean <- taxonomy
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == "") {
    tax.clean[i, 2:7] <- paste0("uncl_", tax.clean[i,1])
  } else if (tax.clean[i,3] == "") {
    tax.clean[i, 3:7] <- paste0("uncl_", tax.clean[i,2])
  } else if (tax.clean[i,4] == "") {
    tax.clean[i, 4:7] <- paste0("uncl_", tax.clean[i,3])
  } else if (tax.clean[i,5] == "") {
    tax.clean[i, 5:7] <- paste0("uncl_", tax.clean[i,4])
  } else if (tax.clean[i,6] == "") {
    tax.clean[i, 6:7] <- paste0("uncl_", tax.clean[i,5])
  } else if (tax.clean[i,7] == "") {
    tax.clean[i,7] <- paste0("uncl_", tax.clean[i,6])
  }
}

# Color palette
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
  '#ffff99','#b15928',"#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD","#AD6F3B",
  "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

# Barplot function
generate_taxa_barplots <- function(my_column, ordered_levels, tax_levels = c("Phylum", "Family", "Genus")) {
  for (ml in tax_levels) {
    taxa.summary <- physeq %>%
      tax_glom(taxrank = ml, NArm = FALSE) %>%
      transform_sample_counts(function(x) x / sum(x)) %>%
      psmelt() %>%
      group_by(.data[[my_column]], .data[[ml]]) %>%
      summarise(Abundance.average = mean(Abundance), .groups = "drop") %>%
      rename(Group = 1, Taxa = 2)
    
    taxa.filtered <- taxa.summary %>%
      group_by(Taxa) %>%
      mutate(overall.max = max(Abundance.average)) %>%
      filter(overall.max > 0.02)
    
    taxa.filtered$Group <- factor(taxa.filtered$Group, levels = ordered_levels)
    taxa.filtered$Taxa <- factor(taxa.filtered$Taxa, levels = names(sort(tapply(taxa.filtered$Abundance.average, taxa.filtered$Taxa, max), decreasing = TRUE)))
    
    p <- ggplot(taxa.filtered, aes(x = Group, y = Abundance.average, fill = Taxa)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = my_colors) +
      ylim(0, 1) +
      labs(y = "Relative Abundance", x = my_column, title = paste0(ml, " (>2%) in at least 1 sample")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
    
    ggsave(paste0("output/taxa/", ml, "_BarPlot_", my_column, ".png"), p, height = 5, width = 5)
  }
}

# Run barplot generation
generate_taxa_barplots("Intervention_ord", c("Baseline", "MedA - Post", "MedWL - Post"))
generate_taxa_barplots("MedDt_ord", c("1. MedDt-WL", "2. MedDt-A"))
generate_taxa_barplots("Obese_ord", c("Overweight", "Obese"))

# --- DESeq2 Analysis ---
# --------------------------------------------
# DESeq2 Analysis (Repaired & Working Version)
# --------------------------------------------

# 1. Set reference levels (IMPORTANT for DESeq2 contrast to work)
metadata$Intervention_ord <- relevel(factor(metadata$Intervention_ord), ref = "Baseline")
metadata$MedDt_ord <- relevel(factor(metadata$MedDt_ord), ref = "1. MedDt-WL")
metadata$Obese_ord <- relevel(factor(metadata$Obese_ord), ref = "Overweight")

# Update metadata in phyloseq object
sample_data(physeq)$Intervention_ord <- metadata$Intervention_ord
sample_data(physeq)$MedDt_ord <- metadata$MedDt_ord
sample_data(physeq)$Obese_ord <- metadata$Obese_ord

# 2. Clean OTU table for DESeq2
OTU.clean <- as.data.frame(otu_table(physeq)) + 1
OTU.physeq <- otu_table(as.matrix(OTU.clean), taxa_are_rows = TRUE)
tax.physeq <- tax_table(as.matrix(tax.clean))
meta.physeq <- sample_data(metadata)
physeq_deseq <- phyloseq(OTU.physeq, tax.physeq, meta.physeq)

# Updated DESeq2 plotting function with top/bottom group labels
run_deseq2 <- function(factor_var, group1, group2, alpha = 0.05) {
  diagdds <- phyloseq_to_deseq2(physeq_deseq, as.formula(paste("~", factor_var)))
  diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")
  
  contrast <- c(factor_var, group1, group2)
  res <- results(diagdds, contrast = contrast, cooksCutoff = FALSE)
  sigtab <- res[which(res$padj < alpha), ]
  sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
  
  sigtab$Genus <- factor(as.character(sigtab$Genus),
                         levels = names(sort(tapply(sigtab$log2FoldChange, sigtab$Genus, max), decreasing = TRUE)))
  
  max_y <- ceiling(max(sigtab$log2FoldChange, na.rm = TRUE))
  min_y <- floor(min(sigtab$log2FoldChange, na.rm = TRUE))
  
  p <- ggplot(sigtab, aes(x = Genus, y = log2FoldChange, color = Phylum)) +
    geom_point(size = 3) +
    ylab(paste0("log2 Fold Change (", group1, " vs ", group2, ")")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          plot.margin = margin(10, 10, 20, 10)) +
    scale_color_manual(values = my_colors[c(2,4,6,8,10,12,14,16,18,20)]) +
    ggtitle(paste0("Differential Abundance: ", group1, " vs ", group2)) +
    annotate("text", x = 1, y = max_y + 0.5, label = paste(group1, "↑"), hjust = 0, fontface = "bold") +
    annotate("text", x = 1, y = min_y - 0.5, label = paste(group2, "↑"), hjust = 0, fontface = "bold")
  
  ggsave(paste0("output/taxa/DESeq2_", factor_var, "_", group1, "_vs_", group2, ".png"), p, width = 10, height = 5)
}

# 🔁 Rerun all comparisons with new labeling
run_deseq2("Intervention_ord", "MedA_Post", "Baseline")
run_deseq2("Intervention_ord", "MedWL_Post", "Baseline")
run_deseq2("MedDt_ord", "2. MedDt-A", "1. MedDt-WL")
run_deseq2("Obese_ord", "Obese", "Overweight")


write_deseq_table <- function(factor_var, group1, group2, alpha = 0.05) {
  diagdds <- phyloseq_to_deseq2(physeq_deseq, as.formula(paste("~", factor_var)))
  diagdds <- DESeq(diagdds)
  res <- results(diagdds, contrast = c(factor_var, group1, group2), cooksCutoff = FALSE)
  sigtab <- res[which(res$padj < alpha), ]
  sigtab <- cbind(as.data.frame(sigtab), as.data.frame(tax_table(physeq_deseq)[rownames(sigtab), ]))
  write.csv(sigtab, file = paste0("output/taxa/DESeq2_", factor_var, "_", group1, "_vs_", group2, ".csv"), row.names = TRUE)
}
