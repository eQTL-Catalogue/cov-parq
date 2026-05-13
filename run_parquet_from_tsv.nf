nextflow.enable.dsl = 2

params.studies_file = params.containsKey('studies_file') ? params.studies_file : null
params.output_base = params.containsKey('output_base') ? params.output_base : null
params.chrom = params.containsKey('chrom') ? params.chrom : null
params.pos_start = params.containsKey('pos_start') ? params.pos_start : null
params.pos_end = params.containsKey('pos_end') ? params.pos_end : null
params.venv_path = params.venv_path ?: "${workflow.projectDir}/.venv"
params.python_bin = params.python_bin ?: "${params.venv_path}/bin/python"
params.script_path = params.containsKey('script_path') ? params.script_path : "${workflow.projectDir}/create_single_parquet_ext.py"

def buildOutputName(String qtlGroup) {
    def outputName = qtlGroup

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
    tag "${study_name}:${qtl_group}"
    publishDir "${params.output_base}/${study_name}", mode: 'copy'

    input:
    tuple val(study_name), val(qtl_group), path(bigwig_dir)

    output:
    path("*.parquet")

    script:
    def outputFileName = buildOutputName(qtl_group)
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
    echo "Study: ${study_name}"
    echo "QTL group: ${qtl_group}"
    echo "Input bigwig dir: ${bigwig_dir}"
    echo "Output file: ${outputFileName}"
    echo "Python: ${params.python_bin}"

    if [ ! -x "${params.python_bin}" ]; then
      echo "ERROR: local venv python not found or not executable at ${params.python_bin}" >&2
      echo "Create it with: python3 -m venv ${params.venv_path} && ${params.venv_path}/bin/pip install -r ${workflow.projectDir}/requirements.txt" >&2
      exit 1
    fi

    "${params.python_bin}" "${params.script_path}" \
      "${bigwig_dir}" \
      "${outputFileName}" \
      --num_processes $task.cpus \
      ${regionArgs.join(' ')}
    """
}

workflow {
    if (!params.studies_file) {
        error "Parameter --studies_file is required."
    }
    if (!params.output_base) {
        error "Parameter --output_base is required."
    }
    if ((params.pos_start != null || params.pos_end != null) && !params.chrom) {
        error "Parameters --pos_start and --pos_end require --chrom."
    }

    def studiesCh = channel
        .fromPath(params.studies_file, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.study_name as String, row.qtl_group as String, file(row.bigwig_dir as String))
        }
        .ifEmpty { error "No rows found in studies file: ${params.studies_file}" }

    RUN_PARQUET(studiesCh)
}
