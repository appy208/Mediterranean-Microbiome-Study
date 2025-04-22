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
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --m-metadata-column Intervention_Status \
  --o-visualization bray-curtis-intervention-significance.qzv \
  --p-method permanova \
  --p-pairwise

echo "🔹 Running PERMANOVA tests for Weighted UniFrac..."

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --m-metadata-column Intervention_Status \
  --o-visualization weighted-unifrac-intervention-significance.qzv \
  --p-method permanova \
  --p-pairwise

echo "🔹 Running PERMANOVA tests for Unweighted UniFrac..."

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --m-metadata-column Intervention_Status \
  --o-visualization unweighted-unifrac-intervention-significance.qzv \
  --p-method permanova \
  --p-pairwise


echo "🔹 Running Alpha Diversity Significance Tests..."

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --m-metadata-column Intervention_Status \
  --o-visualization shannon-intervention-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --m-metadata-column Intervention_Status \
  --o-visualization observed-features-intervention-significance.qzv

echo "🔹 Creating PCoA Emperor Plots..."

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --o-visualization bray-curtis-emperor-intervention.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --o-visualization weighted-unifrac-emperor-intervention.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file Final_Metadata_With_Intervention_Status.txt \
  --o-visualization unweighted-unifrac-emperor-intervention.qzv
