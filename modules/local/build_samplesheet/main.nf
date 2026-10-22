// cohorts -> one sample sheet for this project and sequencing type; see bin/build_samplesheet.py
process BUILD_SAMPLESHEET {
    tag "${project}"

    input:
    val project
    val sequencing_type
    val cohorts

    output:
    path 'samplesheet.tsv'

    script:
    """
    build_samplesheet.py \\
        --project ${project} \\
        --sequencing-type ${sequencing_type} \\
        --cohorts ${cohorts.join(' ')} \\
        --output samplesheet.tsv
    """

    stub:
    """
    printf 'sg_id\\tkind\\tn_units\\tunit\\tr1\\tr2\\treads\\treads_index\\treads_ref\\n' > samplesheet.tsv
    """
}
