// samtools stats of a CRAM, as cpg-flow's CramQcSamtoolsStats
process SAMTOOLS_STATS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(cram), path(crai)
    // the CRAM's own reference, staged apart so it can't collide with the CRAM's name
    tuple path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*'), path(dict, stageAs: 'ref/*')

    output:
    tuple val(meta), path("${meta.id}.samtools-stats")

    script:
    """
    samtools stats \\
        -@${task.cpus - 1} \\
        --reference ${fasta} \\
        ${cram} > ${meta.id}.samtools-stats
    """

    stub:
    """
    touch ${meta.id}.samtools-stats
    """
}
