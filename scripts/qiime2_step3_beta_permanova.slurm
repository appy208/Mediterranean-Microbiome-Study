#!/bin/bash
#SBATCH --account=microbiome
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=qiime2_beta_stats
#SBATCH --output=qiime2_beta_stats_%j.log

echo "🔹 Loading QIIME 2 modules..."
module load biocontainers
module load qiime2

echo "🔹 Changing to working directory..."
cd /scratch/bell/chhetri/ProjectAssignments/data

echo "🔹 Running PERMANOVA tests for Bray-Curtis..."

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column BMI_Category \
  --o-visualization bray-curtis-bmi-significance.qzv \
  --p-method permanova \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column Randomization\ group \
  --o-visualization bray-curtis-randomization-significance.qzv \
  --p-method permanova \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column Time\ Point \
  --o-visualization bray-curtis-timepoint-significance.qzv \
  --p-method permanova \
  --p-pairwise

echo "🔹 Running PERMANOVA tests for Weighted UniFrac..."

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column BMI_Category \
  --o-visualization weighted-unifrac-bmi-significance.qzv \
  --p-method permanova \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column Randomization\ group \
  --o-visualization weighted-unifrac-randomization-significance.qzv \
  --p-method permanova \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column Time\ Point \
  --o-visualization weighted-unifrac-timepoint-significance.qzv \
  --p-method permanova \
  --p-pairwise

echo "🔹 Running PERMANOVA tests for Unweighted UniFrac..."

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column BMI_Category \
  --o-visualization unweighted-unifrac-bmi-significance.qzv \
  --p-method permanova \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column Randomization\ group \
  --o-visualization unweighted-unifrac-randomization-significance.qzv \
  --p-method permanova \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Diet_And_BMI.txt \
  --m-metadata-column Time\ Point \
  --o-visualization unweighted-unifrac-timepoint-significance.qzv \
  --p-method permanova \
  --p-pairwise
