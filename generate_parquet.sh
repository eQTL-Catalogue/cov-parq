#!/bin/bash

#SBATCH --job-name=MacroMap_parq_LCL
#SBATCH --output=MacroMap_parq_LCL.out
#SBATCH --error=MacroMap_parq_LCL.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=01:00:00

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq/MacroMap/LCL/bigwig"
OUTPUT_PATH="/gpfs/helios/projects/eQTLCatalogue/coverage_parquet/MacroMap/LCL_chr22_18500000_20000000.parquet"
PYTHON_BIN=""

echo "Job started on $(hostname) at $(date)"
echo "Study: MacroMap"
echo "Cell type: LCL"
echo "Chromosome filter: 22"
echo "Input: ${INPUT_DIR}"
echo "Output: ${OUTPUT_PATH}"

mkdir -p "$(dirname "${OUTPUT_PATH}")"

PYTHON_BIN="$(command -v python3 || true)"

if [ -z "${PYTHON_BIN}" ] || [ ! -x "${PYTHON_BIN}" ]; then
    echo "Python executable not found. Ensure python3 is on PATH."
    exit 1
fi

echo "Python binary: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c "import sys; print('sys.executable:', sys.executable)"
"${PYTHON_BIN}" -c "import pandas as pd; print('pandas version:', pd.__version__)"

"${PYTHON_BIN}" "${BASE_DIR}/create_single_parquet_ext.py" "${INPUT_DIR}" "${OUTPUT_PATH}" --num_processes 16 --chrom 22 --pos_start 18500000 --pos_end 20000000

echo "Job finished at $(date)"



