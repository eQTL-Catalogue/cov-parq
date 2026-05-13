# cov-parq

**cov-parq** is a Nextflow pipeline designed to efficiently convert large-scale RNA-seq coverage data (in BigWig format) into optimized Parquet files. It processes data per study and subdirectory (e.g., sample groups or conditions), enabling scalable analysis of genomic coverage.

## Features

- **Scalable Processing**: Utilizes Nextflow for parallel execution across samples and subdirectories.
- **Memory Efficient**: Implements a two-stage map-reduce strategy in Python to handle large datasets on systems with limited memory (e.g., HPC nodes).
- **Optimized Output**: Produces sorted Parquet files partitioned by chromosome and genomic position, facilitating fast querying.
- **Flexible Filtering**: Supports filtering by chromosome and specific genomic regions (start/end positions).

## Prerequisites

- **Nextflow**: Ensure Nextflow is installed and available in your path.
- **Python 3**: Python 3.8+ is recommended.
- **Java**: Required for Nextflow (Java 8 or later).

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/cov-parq.git
    cd cov-parq
    ```

2.  **Set up the Python environment:**
    The pipeline expects a local virtual environment named `.venv` in the project root.

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
    ```

## Usage

Run the pipeline using Nextflow. The main entry point is `run_parquet_per_subdir.nf`.

### Basic Command

```bash
nextflow run run_parquet_per_subdir.nf \
    --study <STUDY_NAME> \
    --input_base <PATH_TO_INPUT_DATA> \
    --output_base <PATH_TO_OUTPUT_DIR>
```

For multi-study runs, provide a file instead of `--study`:

```bash
nextflow run run_parquet_per_subdir.nf \
    --studies_file <PATH_TO_STUDY_LIST.tsv> \
    --input_base <PATH_TO_INPUT_DATA> \
    --output_base <PATH_TO_OUTPUT_DIR>
```

### Example

```bash
nextflow run run_parquet_per_subdir.nf \
    --study MacroMap \
    --input_base /gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq \
    --output_base ./output
```

Example multi-study list file (`studies.tsv`, one study per row; first column is used):

```text
MacroMap
GTEx
eQTLGen
```

Example multi-study run:

```bash
nextflow run run_parquet_per_subdir.nf \
    --studies_file studies.tsv \
    --input_base /gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq \
    --output_base ./output
```

### SLURM Execution (HPC)

For running on a SLURM cluster, you can use the provided configuration profile `tartu_hpc` or a submission script like `run_test.sh`.

```bash
# Using the profile directly
nextflow run run_parquet_per_subdir.nf -profile tartu_hpc ...

# Or using a submission script
sbatch run_test.sh
```

## Parameters

| Parameter | Description | Default |
| :--- | :--- | :--- |
| `--study` | Name of a single study to process. Mutually exclusive with `--studies_file`. | N/A |
| `--studies_file` | Path to a text/TSV/CSV file with study names, one per row (first column used). Mutually exclusive with `--study`. Empty lines and `#` comments are ignored. | N/A |
| `--input_base` | Base directory containing study data. Expects structure: `${input_base}/${study}/*/bigwig` | `/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq` |
| `--output_base` | **Required**. Base directory for output Parquet files. | N/A |
| `--chrom` | Optional. Filter by chromosome (e.g., `22`, `X`). | All autosomes (1-22) + X |
| `--pos_start` | Optional. Filter by 1-based inclusive start position (requires `--chrom`). | None |
| `--pos_end` | Optional. Filter by 1-based inclusive end position (requires `--chrom`). | None |
| `--num_processes` | Number of parallel processes per task. | 16 |

## Output

The pipeline generates Parquet files in the specified output directory:
`${output_base}/${study}/<subdir>.parquet`

The output schema includes:
- `seqnames`: Chromosome name (e.g., "1", "X")
- `start`: 1-based start position
- `end`: 1-based end position
- `strand`: Strand information (default `*`)
- `score`: Coverage score
- `sample_id`: Sample identifier derived from the BigWig filename

## Standalone Usage

The core conversion logic is contained in `create_single_parquet_ext.py` and can be run independently:

```bash
python create_single_parquet_ext.py \
    <input_dir> \
    <output_path> \
    --num_processes 8 \
    --chrom 22
```
