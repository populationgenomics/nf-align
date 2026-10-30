// somalier fingerprint of a CRAM, for sample identity and relatedness checks (somalier relate)
process SOMALIER_EXTRACT {
    tag "${meta.id}"

    input:
    // the CRAM decodes only against the reference it was written with: ref_fasta for a new CRAM, or an already
    // aligned CRAM's own. Staged apart, so a FASTA can't collide with the CRAM's name
    tuple val(meta), path(cram), path(crai), path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*')
    path sites // somalier sites VCF

    output:
    tuple val(meta), path("${meta.id}.cram.somalier")

    script:
    """
    somalier extract -d extracted/ --sites ${sites} -f ${fasta} ${cram}
    mv extracted/*.somalier ${meta.id}.cram.somalier
    """

    stub:
    """
    touch ${meta.id}.cram.somalier
    """
}
