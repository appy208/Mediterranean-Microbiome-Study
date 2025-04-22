#!/bin/sh -l
set -e

#SBATCH --nodes=1                     ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 12		              ##number of requested CPUs
#SBATCH --time=04:00:00               ##time requested
#SBATCH --job-name soil_pipeline1     ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using.

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /depot/microbiome/data/ANSC516/ApekshyaChhetri/ProjectAssignments/data/raw

pwd

#step 2, Load Modules

echo "loading biocontainers"
module load biocontainers

echo "loading qiime2-2024.2"
module load qiime2/2024.2


#step 3, Run the code you want

echo "initiating qiime2 moving pictures tutorial pipeline"

echo "downloading necessary files: metadata, sequence files and classifier"

wget \
  -O "sample-metadata.tsv" \
  "https://data.qiime2.org/2023.9/tutorials/atacama-soils/sample_metadata.tsv"

mkdir single-end-sequences

wget \
  -O "emp-paired-end-sequences/forward.fastq.gz" \
  "https://data.qiime2.org/2023.9/tutorials/atacama-soils/10p/forward.fastq.gz"

wget \
  -O "emp-paired-end-sequences/reverse.fastq.gz" \
  "https://data.qiime2.org/2023.9/tutorials/atacama-soils/10p/reverse.fastq.gz"

wget \
  -O "emp-paired-end-sequences/barcodes.fastq.gz" \
  "https://data.qiime2.org/2023.9/tutorials/atacama-soils/10p/barcodes.fastq.gz"


#Get the databases to classify the sequences
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2023.9/common/gg-13-8-99-515-806-nb-classifier.qza"


echo "Data Import"
#Make the sequences into a qiime recognizable format


# Step 1: Import single-end data using a manifest file
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path Manifest-ProjectData.txt \
  --output-path demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Step 2: Summarize imported sequences
echo "Demultiplex."
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# (After this, check quality using demux.qzv and decide trim/trunc values)

# Step 3: Denoise using DADA2 (Example trunc length set to 150; adjust if needed)

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 150 \
  --p-n-threads 0 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats denoising-stats.qza

# Step 4: Summarize feature table
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file Manifest-ProjectData.txt

# Step 5: Visualize representative sequences
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# Step 6: Download pre-trained SILVA
wget https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza
 
# Step 7: Using pre-trained SILVA classifier on qiime to create a taxonomy file
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza


# Step 8: Assign taxonomy 
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza



# END OF SCRIPT




echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch qiime2_step1.slurm
