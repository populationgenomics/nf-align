// record an SG's published CRAM as one Metamist 'cram' analysis, in the SG's own project (meta.project), with its
// index and somalier fingerprint as secondary files and the reference it was written against (meta.ref) in its
// meta. The URIs come from main.nf, built from the same meta fields publishDir publishes with
process REGISTER_ANALYSIS {
    tag "${meta.id}"

    input:
    tuple val(meta), val(cram), val(somalier)
    val sequencing_type
    val ar_guid
    val wait_minutes // how long to wait for the CRAM, index and fingerprint to be published

    output:
    val meta

    script:
    def attempts = wait_minutes * 6
    """
    # publishDir copies asynchronously, so the outputs may not be in GCS yet
    for i in \$(seq ${attempts}); do
        python3 -c "import sys; from cpg_utils.existence_checks import exists; sys.exit(not all(map(exists, sys.argv[1:])))" \\
            ${cram} ${cram}.crai ${somalier} && break
        sleep 10
    done

    # fails (listing them) if any file is still missing
    python3 -m cpg_utils.metamist_registration \\
        --project ${meta.project} \\
        --output ${cram} \\
        --type cram \\
        --secondary index=${cram}.crai somalier=${somalier} \\
        --sgs ${meta.id} \\
        --meta sequencing_type=${sequencing_type} ar-guid=${ar_guid} ref_fasta=${meta.ref}
    """

    stub:
    """
    echo python3 -m cpg_utils.metamist_registration --dry --project ${meta.project} --output ${cram} \\
        --type cram --secondary index=${cram}.crai somalier=${somalier} --sgs ${meta.id} \\
        --meta sequencing_type=${sequencing_type} ar-guid=${ar_guid} ref_fasta=${meta.ref}
    """
}
