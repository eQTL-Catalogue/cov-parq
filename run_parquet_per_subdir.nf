nextflow.enable.dsl = 2

params.input_base = params.input_base ?: '/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq'
params.output_base = params.containsKey('output_base') ? params.output_base : null
params.study = params.containsKey('study') ? params.study : null
params.chrom = params.containsKey('chrom') ? params.chrom : null
params.pos_start = params.containsKey('pos_start') ? params.pos_start : null
params.pos_end = params.containsKey('pos_end') ? params.pos_end : null
params.num_processes = params.num_processes ?: 16
params.venv_path = params.venv_path ?: "${workflow.projectDir}/.venv"
params.python_bin = params.python_bin ?: "${params.venv_path}/bin/python"
params.script_path = params.containsKey('script_path') ? params.script_path : "${workflow.projectDir}/create_single_parquet_ext.py"

def buildOutputName(String subdir) {
    def outputName = subdir

    if (params.chrom) {
        outputName += "_chr${params.chrom}"
    }
    if (params.pos_start != null) {
        outputName += "_${params.pos_start}"
    }
    if (params.pos_end != null) {
        outputName += "_${params.pos_end}"
    }

    return "${outputName}.parquet"
}

process EXTRACT_PARQUET_PARTS {
    tag "$subdir"

    input:
    tuple val(subdir), path(bigwig_dir)

    output:
    tuple val(subdir), path("${subdir}_parquet_parts")

    script:
    def regionArgs = []
    if (params.chrom) {
        regionArgs << "--chrom ${params.chrom}"
    }
    if (params.pos_start != null) {
        regionArgs << "--pos_start ${params.pos_start}"
    }
    if (params.pos_end != null) {
        regionArgs << "--pos_end ${params.pos_end}"
    }

    """
    echo "Extracting intermediate Parquet parts for subdir: ${subdir}"
    echo "Input: ${bigwig_dir}"
    echo "Parts directory: ${subdir}_parquet_parts"
    echo "Python: ${params.python_bin}"

    if [ ! -x "${params.python_bin}" ]; then
      echo "ERROR: local venv python not found or not executable at ${params.python_bin}" >&2
      echo "Create it with: python3 -m venv ${params.venv_path} && ${params.venv_path}/bin/pip install -r ${workflow.projectDir}/requirements.txt" >&2
      exit 1
    fi

    "${params.python_bin}" "${params.script_path}" \
      "${bigwig_dir}" \
      "${subdir}_parquet_parts" \
      --mode extract \
      --num_processes ${params.num_processes} \
      ${regionArgs.join(' ')}
    """
}

process MERGE_PARQUET_PARTS {
    tag "$subdir"
    publishDir "${params.output_base}/${params.study}", mode: 'copy'

    input:
    tuple val(subdir), path(parts_dir)

    output:
    path("*.parquet")

    script:
    def outputFileName = buildOutputName(subdir)
    def publishDestination = "${params.output_base}/${params.study}"

    def regionArgs = []
    if (params.chrom) {
        regionArgs << "--chrom ${params.chrom}"
    }
    if (params.pos_start != null) {
        regionArgs << "--pos_start ${params.pos_start}"
    }
    if (params.pos_end != null) {
        regionArgs << "--pos_end ${params.pos_end}"
    }

    """
    echo "Merging Parquet parts for subdir: ${subdir}"
    echo "Parts directory: ${parts_dir}"
    echo "Local output file: ${outputFileName}"
    echo "Publish destination: ${publishDestination}"
    echo "Python: ${params.python_bin}"

    if [ ! -x "${params.python_bin}" ]; then
      echo "ERROR: local venv python not found or not executable at ${params.python_bin}" >&2
      echo "Create it with: python3 -m venv ${params.venv_path} && ${params.venv_path}/bin/pip install -r ${workflow.projectDir}/requirements.txt" >&2
      exit 1
    fi

    "${params.python_bin}" "${params.script_path}" \
      "${parts_dir}" \
      "${outputFileName}" \
      --mode merge \
      --num_processes ${params.num_processes} \
      ${regionArgs.join(' ')}
    """
}

workflow {
    if (!params.output_base) {
        error "Parameter --output_base is required."
    }

    if (!params.study) {
        error "Parameter --study is required."
    }

    if ((params.pos_start != null || params.pos_end != null) && !params.chrom) {
        error "Parameters --pos_start and --pos_end require --chrom."
    }

    def inputPattern = "${params.input_base}/${params.study}/*/bigwig"

    channel
        .fromPath(inputPattern, type: 'dir', checkIfExists: true)
        .map { bigwigDir -> tuple(bigwigDir.parent.name, bigwigDir) }
        .map { subdir, _bigwigDir -> "Will write: ${params.output_base}/${params.study}/${buildOutputName(subdir)}" }
        .toList()
        .view { lines ->
            "Found ${lines.size()} subdirectories with bigwig data.\n${lines.join('\n')}"
        }

    def bigwigDirsCh = channel
        .fromPath(inputPattern, type: 'dir', checkIfExists: true)
        .ifEmpty { error "No input directories found for pattern: ${inputPattern}" }
        .map { bigwigDir -> tuple(bigwigDir.parent.name, bigwigDir) }

    def parquetPartsCh = EXTRACT_PARQUET_PARTS(bigwigDirsCh)
    MERGE_PARQUET_PARTS(parquetPartsCh)
}
