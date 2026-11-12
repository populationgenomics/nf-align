/*
 * Record one published file as a Metamist analysis of `type`, in the SG's own project (meta.project), with its
 * secondary files and analysis meta. Included as REGISTER_CRAM (the CRAM) and REGISTER_CRAM_QC (its
 * QC metrics). The URIs come from main.nf, built from the same meta fields publishDir publishes with
 */
process REGISTER_ANALYSIS {
    tag "${meta.id}"

    input:
    // secondary: [key: uri]; analysis_meta: [key: value]. Both keep their order on the command line
    tuple val(meta), val(output), val(secondary), val(analysis_meta)
    val type
    val wait_minutes // how long to wait for every file to be published

    output:
    val meta

    script:
    def attempts = wait_minutes * 6
    def files = ([output] + secondary.values()).join(' ')
    def args = registrationArgs(meta, output, type, secondary, analysis_meta)
    """
    # publishDir copies asynchronously, so the outputs may not be in GCS yet
    for i in \$(seq ${attempts}); do
        python3 -c "import sys; from cpg_utils.existence_checks import exists; sys.exit(not all(map(exists, sys.argv[1:])))" \\
            ${files} && break
        sleep 10
    done

    # fails (listing them) if any file is still missing
    python3 -m cpg_utils.metamist_registration ${args}
    """

    stub:
    def args = registrationArgs(meta, output, type, secondary, analysis_meta)
    """
    echo python3 -m cpg_utils.metamist_registration --dry ${args}
    """
}

def registrationArgs(meta, output, type, secondary, analysis_meta) {
    "--project ${meta.project} --output ${output} --type ${type} --secondary ${keyValues(secondary)} --sgs ${meta.id} --meta ${keyValues(analysis_meta)}"
}

def keyValues(map) {
    map.collect { key, value -> "${key}=${value}" }.join(' ')
}
