// record a published CRAM as a Metamist 'cram' analysis, in the SG's own project (meta.project). cram_uri
// comes from destination() in main.nf, built from the same meta fields publishDir publishes with
process REGISTER_ANALYSIS {
    tag "${meta.id}"

    input:
    tuple val(meta), val(cram_uri)
    val sequencing_type
    val ar_guid
    val fasta
    val wait_minutes // how long to wait for the CRAM and index to be published

    output:
    val meta

    script:
    def attempts = wait_minutes * 6
    """
    # publishDir copies asynchronously, so the CRAM and index may not be in GCS yet
    for i in \$(seq ${attempts}); do
        python3 -c "import sys; from cpg_utils.existence_checks import exists; sys.exit(not all(map(exists, sys.argv[1:])))" \\
            ${cram_uri} ${cram_uri}.crai && break
        sleep 10
    done

    # fails (listing them) if either file is still missing
    python3 -m cpg_utils.metamist_registration \\
        --project ${meta.project} \\
        --output ${cram_uri} \\
        --type cram \\
        --secondary index=${cram_uri}.crai \\
        --sgs ${meta.id} \\
        --meta sequencing_type=${sequencing_type} ar-guid=${ar_guid} ref_fasta=${fasta}
    """

    stub:
    """
    echo python3 -m cpg_utils.metamist_registration --dry --project ${meta.project} --output ${cram_uri} \\
        --type cram --secondary index=${cram_uri}.crai --sgs ${meta.id} \\
        --meta sequencing_type=${sequencing_type} ar-guid=${ar_guid} reference=${fasta}
    """
}
