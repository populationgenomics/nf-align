// record a published CRAM as a Metamist 'cram' analysis; cram_uri comes from cramPath() in main.nf,
// the same function the output block publishes with
process REGISTER_ANALYSIS {
    tag "${meta.id}"

    input:
    tuple val(meta), val(cram_uri)
    val project

    output:
    val meta

    script:
    def attempts = params.register_wait_minutes.toInteger() * 6
    """
    # workflow outputs are published asynchronously, so the CRAM and index may not be in GCS yet
    for i in \$(seq ${attempts}); do
        python3 -c "import sys; from cpg_utils.existence_checks import exists; sys.exit(not all(map(exists, sys.argv[1:])))" \\
            ${cram_uri} ${cram_uri}.crai && break
        sleep 10
    done

    # fails (listing them) if either file is still missing
    python3 -m cpg_utils.metamist_registration \\
        --project ${project} \\
        --output ${cram_uri} \\
        --type cram \\
        --secondary index=${cram_uri}.crai \\
        --sgs ${meta.id} \\
        --meta sequencing_type=${params.sequencing_type} ar-guid=${params.arGuid}
    """

    stub:
    """
    echo python3 -m cpg_utils.metamist_registration --dry --project ${project} --output ${cram_uri} \\
        --type cram --secondary index=${cram_uri}.crai --sgs ${meta.id} \\
        --meta sequencing_type=${params.sequencing_type} ar-guid=${params.arGuid}
    """
}
