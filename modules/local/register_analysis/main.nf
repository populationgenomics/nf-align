// record a published output as a Metamist analysis, in the SG's own project (meta.project): a 'cram' with its
// index, or a 'somalier' fingerprint. output_uri comes from destination() in main.nf, built from the same meta
// fields publishDir publishes with
process REGISTER_ANALYSIS {
    tag "${meta.id}:${type}"

    input:
    tuple val(meta), val(type), val(output_uri), val(fasta) // fasta: the reference the CRAM was written against
    val sequencing_type
    val ar_guid
    val wait_minutes // how long to wait for the output (and index) to be published

    output:
    tuple val(meta), val(type)

    script:
    def attempts = wait_minutes * 6
    // a CRAM is registered with its index, which has to be published too
    def files = type == 'cram' ? "${output_uri} ${output_uri}.crai" : output_uri
    def secondary = type == 'cram' ? "--secondary index=${output_uri}.crai" : ''
    """
    # publishDir copies asynchronously, so the output may not be in GCS yet
    for i in \$(seq ${attempts}); do
        python3 -c "import sys; from cpg_utils.existence_checks import exists; sys.exit(not all(map(exists, sys.argv[1:])))" \\
            ${files} && break
        sleep 10
    done

    # fails (listing them) if any file is still missing
    python3 -m cpg_utils.metamist_registration \\
        --project ${meta.project} \\
        --output ${output_uri} \\
        --type ${type} \\
        ${secondary} \\
        --sgs ${meta.id} \\
        --meta sequencing_type=${sequencing_type} ar-guid=${ar_guid} ref_fasta=${fasta}
    """

    stub:
    def secondary = type == 'cram' ? "--secondary index=${output_uri}.crai " : ''
    """
    echo python3 -m cpg_utils.metamist_registration --dry --project ${meta.project} --output ${output_uri} \\
        --type ${type} ${secondary}--sgs ${meta.id} \\
        --meta sequencing_type=${sequencing_type} ar-guid=${ar_guid} ref_fasta=${fasta}
    """
}
