#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="create_cov_parquet"

module load any/jdk/1.8.0_265
module load nextflow
module load squashfs/4.4

nextflow run run_parquet_per_subdir.nf -profile tartu_hpc -resume \
  --study MacroMap \
  --input_base /gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq \
  --output_base /gpfs/helios/projects/eQTLCatalogue/coverage_parquet/MacroMap