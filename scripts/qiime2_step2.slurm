#!/bin/sh -l
set -e  # Exit script if any command fails

#SBATCH --nodes=1                 ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 12		           ##number of requested CPUs
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name soil_pipeline2    ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using.

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /depot/microbiome/data/ANSC516/ApekshyaChhetri/ProjectAssignments/data/output

pwd

#step 2, Load Modules

echo "loading biocontainers"
module load biocontainers

echo "loading qiime2-2024.2"
module load qiime2/2024.2

#step 3, Run the code you want

qiime tools export \
--input-path demux-subsample.qzv \
--output-path ./demux-subsample/

qiime demux filter-samples \
--i-demux demux-subsample.qza \
--m-metadata-file ./demux-subsample/per-sample-fastq-counts.tsv \
--p-where 'CAST([forward sequence count] AS INT) > 100' \
--o-filtered-demux demux.qza


##De-noise. Error corrections.

echo "dada2"

qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux.qza \
--p-trim-left-f 13 \
--p-trim-left-r 13 \
--p-trunc-len-f 150 \
--p-trunc-len-r 150 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza

qiime feature-classifier classify-sklearn \
--i-classifier gg-13-8-99-515-806-nb-classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

#Convert .qza to .qzv format

qiime metadata tabulate \
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv


#Create a .qzv from the output of denoise
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv


echo "Create a phylogenetic tree. (Phylogenetic method)"
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

echo "Rarefaction"
qiime diversity alpha-rarefaction \
--i-table table.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 2000 \
--m-metadata-file sample-metadata.tsv \
--o-visualization alpha-rarefaction.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch qiime2_step2.slurm
