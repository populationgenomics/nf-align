// cohorts or SG IDs, from any datasets at this access level -> one sample sheet for this sequencing type;
// see bin/build_samplesheet.py
process BUILD_SAMPLESHEET {
    tag "${access_level}"

    input:
    val access_level
    val sequencing_type
    val cohorts
    val sgids
    val run_dir // only read by publishDir (conf/base.config)

    output:
    path 'samplesheet.tsv'

    script:
    // exactly one of cohorts and sgids is non-empty
    def ids = cohorts ? "--cohorts ${cohorts.join(' ')}" : "--sgids ${sgids.join(' ')}"
    """
    build_samplesheet.py \\
        --access-level ${access_level} \\
        --sequencing-type ${sequencing_type} \\
        ${ids} \\
        --output samplesheet.tsv
    """

    stub:
    """
    printf 'sg_id\\tdataset\\tkind\\tn_units\\tunit\\tr1\\tr2\\treads\\treads_index\\treads_ref\\n' > samplesheet.tsv
    """
}
