# cov-parq

`cov-parq` converts RNA-seq coverage files from BigWig (`.bw` or `.bigwig`) into sorted Parquet files. The conversion is designed for large eQTL Catalogue-style study folders and can be run either through Nextflow on SLURM or directly with Python for a single BigWig directory.

The Python converter uses a two-stage map-reduce approach:

1. Each BigWig file is converted to an intermediate Parquet part in parallel.
2. DuckDB merges the parts into one final Parquet file, sorted by chromosome and start position.

## Repository Layout

- `run_parquet_from_tsv.nf`: production Nextflow entry point. Reads a TSV with explicit `study_name`, `qtl_group`, and `bigwig_dir` columns.
- `run_parquet_per_subdir.nf`: convenience Nextflow entry point for one study under `${input_base}/${study}/*/bigwig`.
- `create_single_parquet_ext.py`: standalone BigWig directory to Parquet converter.
- `scripts/find_bigwig_dirs.py`: scans an RNA-seq root and writes a TSV for `run_parquet_from_tsv.nf`.
- `scripts/check_bigwig.py`: checks BigWig readability and prints example intervals.
- `scripts/test_output.py`: validates output Parquet schema, coordinate filters, and row order.
- `scripts/check_single_parquet.py`: checks one Parquet file for null values by row group.

## Requirements

- Java and Nextflow for pipeline runs.
- Python 3 with packages from `requirements.txt`: `duckdb`, `pandas`, `pyarrow`, `pyBigWig`, and `tqdm`.
- On Tartu HPC, the provided run scripts load `any/jdk/1.8.0_265`, `nextflow`, and `squashfs/4.4`.

Create the local Python environment expected by the Nextflow workflows:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Input Layout

The expected source data layout is:

```text
<input_base>/<study_name>/<qtl_group>/bigwig/*.bw
<input_base>/<study_name>/<qtl_group>/bigwig/*.bigwig
```

For TSV-driven runs, the TSV must have a header with these columns:

```text
study_name	qtl_group	bigwig_dir
ROSMAP	brain_naive	/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq/ROSMAP/brain_naive/bigwig
Walker_2019	Neocortex	/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq/Walker_2019/Neocortex/bigwig
```

Generate this TSV from an input root with:

```bash
python scripts/find_bigwig_dirs.py \
  --root /gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq \
  --out studies.auto.tsv
```

## Running With Nextflow

Use `run_parquet_from_tsv.nf` when you already have an explicit study list:

```bash
nextflow run run_parquet_from_tsv.nf -profile tartu_hpc -resume \
  --studies_file studies.tsv \
  --output_base /gpfs/helios/projects/eQTLCatalogue/coverage_parquet
```

The repository also includes `run_multi_study.sh`, which submits the same TSV-driven workflow to SLURM.

Use `run_parquet_per_subdir.nf` for one study when the input follows the standard directory layout:

```bash
nextflow run run_parquet_per_subdir.nf -profile tartu_hpc -resume \
  --study MacroMap \
  --input_base /gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq \
  --output_base /gpfs/helios/home/kerimov/alasoo_lab/cov-parq/test
```

## Region Filtering

By default, the converter keeps autosomes `1` through `22`. Use `--chrom` to keep one chromosome; `X` is supported only when requested explicitly.

Position filters are 1-based inclusive and require `--chrom`:

```bash
nextflow run run_parquet_from_tsv.nf -profile tartu_hpc -resume \
  --studies_file studies.tsv \
  --output_base /path/to/output \
  --chrom 22 \
  --pos_start 18500000 \
  --pos_end 20000000
```

Output filenames include the filter suffix, for example `CIL_6_chr22_18500000_20000000.parquet`.

## Parameters

| Parameter | Used by | Description | Default |
| --- | --- | --- | --- |
| `--studies_file` | `run_parquet_from_tsv.nf` | TSV with `study_name`, `qtl_group`, and `bigwig_dir` columns. | Required |
| `--study` | `run_parquet_per_subdir.nf` | Single study name to scan under `input_base`. | Required |
| `--input_base` | `run_parquet_per_subdir.nf` | Root containing `<study>/<qtl_group>/bigwig` directories. | `/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq` |
| `--output_base` | both workflows | Base directory for published Parquet files. | Required |
| `--chrom` | both workflows | Optional chromosome to keep, e.g. `22`, `chr22`, or `X`. | Autosomes `1-22` |
| `--pos_start` | both workflows | Optional 1-based inclusive start filter. Requires `--chrom`. | None |
| `--pos_end` | both workflows | Optional 1-based inclusive end filter. Requires `--chrom`. | None |
| `--venv_path` | both workflows | Python virtual environment path. | `${workflow.projectDir}/.venv` |
| `--python_bin` | both workflows | Python executable used by Nextflow tasks. | `${venv_path}/bin/python` |
| `--script_path` | both workflows | Converter script path. | `${workflow.projectDir}/create_single_parquet_ext.py` |
| `--num_processes` | `run_parquet_per_subdir.nf` | Parallel workers passed to the converter. | `16` |

For `run_parquet_from_tsv.nf`, the converter uses the Nextflow task CPU count (`$task.cpus`) as `--num_processes`.

## Output

Each input BigWig directory produces one Parquet file:

```text
<output_base>/<study_name>/<qtl_group>.parquet
```

The output schema is:

| Column | Type | Description |
| --- | --- | --- |
| `seqnames` | string | Chromosome without `chr`, for example `1`, `22`, or `X`. |
| `start` | int32 | 1-based inclusive start coordinate. |
| `end` | int32 | 1-based inclusive end coordinate. |
| `strand` | string | Always `*`. |
| `score` | float64 | Coverage score from the BigWig interval. |
| `sample_id` | string | Sample identifier derived from the BigWig filename before the first dot. |

BigWig intervals are converted from 0-based half-open coordinates to 1-based inclusive coordinates.

## Validation And Inspection

Check BigWig readability:

```bash
python scripts/check_bigwig.py /path/to/bigwig_dir
python scripts/check_bigwig.py /path/to/sample.bigwig --chrom chr22 --n 20
```

Validate an output Parquet file:

```bash
python scripts/test_output.py /path/to/CIL_6.parquet
python scripts/test_output.py /path/to/CIL_6_chr22_18500000_20000000.parquet \
  --expected-chrom 22 \
  --expected-pos-start 18500000 \
  --expected-pos-end 20000000
```

Check for null values in a Parquet file:

```bash
python scripts/check_single_parquet.py /path/to/CIL_6.parquet
```

## Standalone Converter

Run the converter directly for one BigWig directory:

```bash
python create_single_parquet_ext.py \
  /path/to/bigwig \
  /path/to/output.parquet \
  --num_processes 8 \
  --chrom 22 \
  --pos_start 18500000 \
  --pos_end 20000000
```
