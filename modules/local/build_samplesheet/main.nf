// cohorts or SG IDs -> one sample sheet for this project and sequencing type; see bin/build_samplesheet.py
process BUILD_SAMPLESHEET {
    tag "${project}"

    input:
    val project
    val sequencing_type
    val cohorts
    val sgids

    output:
    path 'samplesheet.tsv'

    script:
    // exactly one of cohorts and sgids is non-empty
    def ids = cohorts ? "--cohorts ${cohorts.join(' ')}" : "--sgids ${sgids.join(' ')}"
    """
    build_samplesheet.py \\
        --project ${project} \\
        --sequencing-type ${sequencing_type} \\
        ${ids} \\
        --output samplesheet.tsv
    """

    stub:
    """
    printf 'sg_id\\tkind\\tn_units\\tunit\\tr1\\tr2\\treads\\treads_index\\treads_ref\\n' > samplesheet.tsv
    """
}
