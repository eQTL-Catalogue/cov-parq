nextflow.enable.dsl = 2

params.input_base = params.input_base ?: '/gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq'
params.output_base = params.output_base ?: '/gpfs/helios/projects/eQTLCatalogue/coverage_parquet'
params.study = params.study ?: 'MacroMap'
params.chrom = params.chrom ?: '22'
params.pos_start = params.pos_start ?: 18500000
params.pos_end = params.pos_end ?: 20000000
params.num_processes = params.num_processes ?: 16
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

process RUN_PARQUET {
    tag "$subdir"

    input:
    tuple val(subdir), path(bigwig_dir)

    script:
    def outputFileName = buildOutputName(subdir)
    def outputPath = "${params.output_base}/${params.study}/${outputFileName}"

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
    echo "Processing subdir: ${subdir}"
    echo "Input: ${bigwig_dir}"
    echo "Output: ${outputPath}"

    python3 "${params.script_path}" \
      "${bigwig_dir}" \
      "${outputPath}" \
      --num_processes ${params.num_processes} \
      ${regionArgs.join(' ')}
    """
}

workflow {
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

    RUN_PARQUET(bigwigDirsCh)
}
