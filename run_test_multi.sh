#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="create_cov_parquet"

module load any/jdk/1.8.0_265
module load nextflow
module load squashfs/4.4

# nextflow run run_parquet_per_subdir.nf -profile tartu_hpc -resume \
#   --study MacroMap \
#   --input_base /gpfs/helios/home/kerimov/alasoo_lab/cov-parq \
#   --output_base /gpfs/helios/home/kerimov/alasoo_lab/cov-parq/test


nextflow run run_parquet_from_tsv.nf -profile tartu_hpc -resume \
  --studies_file studies.tsv \
  --input_base /gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq \
  --output_base /gpfs/helios/home/kerimov/alasoo_lab/cov-parq/test_results