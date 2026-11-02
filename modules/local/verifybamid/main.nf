// VerifyBamID2 contamination estimate (FREEMIX) of a CRAM, as cpg-flow's CramQcVerifyBamId
process VERIFYBAMID {
    tag "${meta.id}"

    input:
    // the CRAM's own reference, staged apart so it can't collide with the CRAM's name
    tuple val(meta), path(cram), path(crai), path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*'), path(dict, stageAs: 'ref/*')
    // the 1000G contamination resources for the run's sequencing type
    tuple path(ud), path(bed), path(mu)
    val sequencing_type

    output:
    tuple val(meta), path("${meta.id}.verify-bamid.selfSM")

    script:
    // cpg-flow caps exome depth, as targeted regions can be very deep
    def max_depth = sequencing_type == 'exome' ? '--max-depth 1000' : ''
    """
    VerifyBamID \\
        --NumThread ${task.cpus} \\
        --Verbose \\
        --NumPC 4 \\
        --Output OUTPUT \\
        --BamFile ${cram} \\
        --Reference ${fasta} \\
        --UDPath ${ud} \\
        --MeanPath ${mu} \\
        --BedPath ${bed} \\
        ${max_depth} 1>/dev/null
    mv OUTPUT.selfSM ${meta.id}.verify-bamid.selfSM
    """

    stub:
    """
    touch ${meta.id}.verify-bamid.selfSM
    """
}
