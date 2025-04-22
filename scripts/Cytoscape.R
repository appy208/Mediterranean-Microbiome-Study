# Genus-Level Co-occurrence Network Graph (Spearman Correlation)
# Author: Apekshya's Project | Created: 2025

# ------------------------------
# STEP 1: Load libraries
# ------------------------------
library(openxlsx)
library(igraph)
library(ggraph)
library(tidyverse)
library(Hmisc)
library(reshape2)
library(qiime2R)



# ------------------------------
# STEP 5: Save CSVs for Cytoscape
# ------------------------------
if (!dir.exists("output/cooccurrence")) dir.create("output/cooccurrence", recursive = TRUE)
write.csv(strong_results_taxa, "output/cooccurrence/strong_results_taxa_all.csv", row.names = FALSE)

for (grp in treatments) {
  df_group <- strong_results_taxa %>%
    filter(trt == grp, abs(rho) >= 0.6, qval < 0.05, Genus.x != Genus.y) %>%
    select(Source = Genus.x, Target = Genus.y, rho, qval, trt)
  write.csv(df_group, paste0("output/cooccurrence/edges_", gsub(" ", "_", grp), ".csv"), row.names = FALSE)
  
  node_list <- unique(c(df_group$Source, df_group$Target))
  write.csv(data.frame(ID = node_list), paste0("output/cooccurrence/nodes_", gsub(" ", "_", grp), ".csv"), row.names = FALSE)
}

# ------------------------------
# STEP 6: Optional R network plot (example: Obese)
# ------------------------------
obese_edges <- strong_results_taxa %>%
  filter(trt == "Obese", abs(rho) >= 0.6, qval < 0.05, Genus.x != Genus.y) %>%
  select(Source = Genus.x, Target = Genus.y, rho)

g_obese <- graph_from_data_frame(obese_edges, directed = FALSE)

graph_plot <- ggraph(g_obese, layout = "fr") +
  geom_edge_link(aes(edge_alpha = abs(rho), edge_width = abs(rho), edge_color = rho)) +
  geom_node_point(size = 5, color = "steelblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_color_gradient2(low = "red", mid = "gray", high = "blue", midpoint = 0) +
  theme_void() +
  labs(title = "Obese: Genus Co-occurrence Network (|rho| ≥ 0.6, q < 0.05)",
       edge_color = "Spearman's rho")

print(graph_plot)
ggsave("output/cooccurrence/genus_network_obese.pdf", plot = graph_plot, width = 10, height = 8)
