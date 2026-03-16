#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --job-name="create_cov_parquet"

module load any/jdk/1.8.0_265
module load nextflow
module load py-pyarrow/10.0.1
module load squashfs/4.4

nextflow run run_parquet_per_subdir.nf -profile tartu_hpc -resume \
  --study MacroMap \
  --input_base /gpfs/helios/home/kerimov/alasoo_lab/cov-parq
