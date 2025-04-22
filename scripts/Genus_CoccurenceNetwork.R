# Genus-Level Co-occurrence Network Graph (Spearman Correlation)
# ------------------------------
# STEP 1: Load libraries
# ------------------------------
library(openxlsx)
library(igraph)
install.packages("ggraph", dependencies = TRUE)
library(ggraph)
library(tidyverse)
library(Hmisc)
library(reshape2)
library(qiime2R)

# ------------------------------
# STEP 2: Load ASV table and metadata
# ------------------------------
# Replace with your own ASV table and metadata paths
setwd("~/Downloads/ANSC 516/Final Project_Apekshya/data")

#load metadata
metadata <- read.delim("Final_Metadata_With_Intervention_Status.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#manual change for data type
metadata$Intervention_Status <- factor(metadata$Intervention_Status, levels = c("Baseline", "MedA - Post", "MedWL - Post"))

#load ASV
ASVs <- read_qza("core-metrics-results/rarefied_table.qza")
ASV_table <- as.data.frame(ASVs$data)

# Add ASV labels
##I added in this chuck. What does it do and why?

#creating a column named ASVnos which assigns the values as: ASVn where n changes with row number
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))

##creating a column named ASVstring which assigns the values as: ASVn where n changes with row number
ASV_table$ASVstring <- rownames(ASV_table)

rownames(ASV_table) <- ASV_table$ASVnos

#creates a new table ASVkey after copying ASVnos and ASVstring
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)]
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
######################################################################

#transforming dataset table by changing rows and columns since R needs ASVs in row for analysis
dataset <- as.data.frame(t(ASV_table))

# we are going to create a network per treatment
head(dataset[,1:10])

# Create grouping variable
##Creating grouping and ordinance based on diet intervention
metadata$Intervention_ord <- factor(metadata$Intervention_Status, 
                                    levels = c("Baseline", "MedA - Post", "MedWL - Post"))

my_column <- "Intervention_ord"
treatments <- unique(dataset[[my_column]])
num_metadata_columns <- ncol(meta_raw)

# Load and join metadata
meta_raw <- read.delim("Final_Metadata_With_Intervention_Status.txt", 
                       header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)

str(meta_raw)
colnames(meta_raw)[4] = "BMI_Category"

# Assign correct column names manually
colnames(meta_raw) <- c("SampleID", "Sample_Name", "BMI", "BMI_Category", "Time_Point", 
                        "Randomization_Group", "Percent_Weight_Change", 
                        "MedDiet_Adherence_Change", "MediScore_Change",
                        "Classes_Attended", "Class_Attendance_Pct", "Age", 
                        "Intervention_Status")

dim(dataset)
head(dataset[, 1:5])


# Sample number thresholds per group (for correlation robustness)
n1 <- 8
n2 <- 8
n3 <- 8
q_cutoff <- 0.05

for(i in 1:length(treatments)){
  #subset the data for a particular treatment YOU MUST ENTER THE HEADER OF THE COLUMN THAT HAS THE DIFFERENT TREATMENTS IN THIS CASE “Foaming_Status”
  print(paste("reading ",treatments[i],sep=""))
  temp<-subset(meta_raw, get(my_column)==treatments[i])
  tempn<-subset(meta_raw, get(my_column)==treatments[i])
  print(paste("finished reading ",treatments[i],sep=""))
  # making an object that has all the results in it (both rho and P values)
  results<-rcorr(as.matrix(temp[,-c(1:num_metadata_columns)]),type="spearman") ## use the "-c" parameter to remove metadata columns
  resultsn<-rcorr(as.matrix(tempn[,-c(1:num_metadata_columns)]),type="spearman")
  
  #make two seperate objects for p-value and correlation coefficients
  rhos<-results$r
  ps<-results$P
  ns<-resultsn$n
  # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
  ps_melt<-na.omit(melt(ps))
  #creating a qvalue based on FDR
  ps_melt$qval<-p.adjust(ps_melt$value, method = "BH")
  #making column names more relevant
  
  names(ps_melt)[3]<-"pval"
  # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
  ps_sub<-subset(ps_melt, qval < q_cutoff)
  
  # now melting the rhos, note the similarity between ps_melt and rhos_melt
  rhos_melt<-na.omit(melt(rhos))
  names(rhos_melt)[3]<-"rho"
  
  # now melting the ns
  ns_melt<-(melt(ns))
  names(ns_melt)[3]<-"n"
  
  #merging together and remove negative rhos
  merged<-merge(ps_sub,rhos_melt,by=c("Var1","Var2"))
  if (treatments[i]==treatments[1]) {
    merged<-merge(merged,subset(ns_melt, n > n1),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[2]) {
    merged<-merge(merged,subset(ns_melt, n > n2),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[3]) {
    merged<-merge(merged,subset(ns_melt, n > n3),by=c("Var1","Var2"))
  }   else
    print("Somethings wrong with your treatment designations. Please Check!!")
  
  if (nrow(merged) > 0) {
    merged$trt<-treatments[i]
    final_results<-rbind(final_results, merged)
  }   else {
    print("no correlations for this variable")
  }
  
  print(paste("finished ",treatments[i],sep=""))
}


##Creating grouping and ordinance based on BMI/Weight
metadata$Weight_ord <- factor(metadata$BMI_Category, 
                              levels = c("Obese","Overweight"))
my_column <- "Weight_ord"
treatments <- unique(dataset[[my_column]])
num_metadata_columns <- ncol(metadata)

# ------------------------------
# STEP 4: Merge taxonomy for Genus-level graph
# ------------------------------
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

# Merge taxonomy with correlation results
strong_results <- final_results
strong_results_taxa <- merge(strong_results, ASVkey, by.x = "Var1", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, ASVkey, by.x = "Var2", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.x", by.y = 0)
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.y", by.y = 0)

# ------------------------------
# STEP 5: Filter and save genus-level correlations
# ------------------------------
df <- strong_results_taxa

df_filtered <- df %>%
  filter(abs(rho) >= 0.6 & qval < 0.05) %>%
  filter(Genus.x != Genus.y) %>%
  select(Genus.x, Genus.y, rho, qval, trt)

write.csv(df_filtered, "output/cooccurrence/filtered_genus_cooccurrence.csv", row.names = FALSE)
write.xlsx(df_filtered, "output/cooccurrence/filtered_genus_cooccurrence.xlsx")

# ------------------------------
# STEP 6: Create network graph
# ------------------------------
g <- graph_from_data_frame(df_filtered[, 1:3], directed = FALSE)

# ------------------------------
# STEP 7: Plot network
# ------------------------------
graph_plot <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_alpha = abs(rho), edge_width = abs(rho), edge_color = rho)) +
  geom_node_point(size = 5, color = "steelblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_color_gradient2(low = "red", mid = "gray", high = "blue", midpoint = 0) +
  theme_void() +
  labs(title = "Genus Co-occurrence Network (|rho| ≥ 0.6, q < 0.05)",
       edge_color = "Spearman's rho")

print(graph_plot)

# ------------------------------
