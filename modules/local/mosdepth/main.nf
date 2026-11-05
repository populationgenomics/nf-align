/*
 * mosdepth coverage of a CRAM: the whole-genome depth distribution, the mean depth in 100 bp windows, and the
 * genome split into depth bands (0, 1-4, 5-19, 20-49, 50+). No per-base output.
 *
 * mosdepth walks the CRAM chromosome by chromosome through its index, so it can't read the decoded stream
 * PICARD_METRICS tees to its JVMs: it decodes the CRAM itself, in its own task.
 */
process MOSDEPTH {
    tag "${meta.id}"

    input:
    // the CRAM's own reference, staged apart so it can't collide with the CRAM's name
    tuple val(meta), path(cram), path(crai), path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*'), path(dict, stageAs: 'ref/*')

    output:
    tuple val(meta),
        path("${meta.id}.mosdepth.global.dist.txt"),
        path("${meta.id}.mosdepth.region.dist.txt"),
        path("${meta.id}.mosdepth.summary.txt"),
        path("${meta.id}.regions.bed.gz"),
        path("${meta.id}.regions.bed.gz.csi"),
        path("${meta.id}.quantized.bed.gz"),
        path("${meta.id}.quantized.bed.gz.csi")

    script:
    // decompression threads, on top of mosdepth's own; the mosdepth docs see no gain past 4
    def threads = Math.max(1, Math.min(4, task.cpus - 1))
    // --fast-mode ignores CIGAR operations and doesn't correct mate overlaps, which makes decoding a CRAM much
    // faster. The cost: a read covers everything from its start to its end, deletions included, and where a
    // pair's mates overlap, each base is counted twice
    """
    mosdepth \\
        --threads ${threads} \\
        --fasta ${fasta} \\
        --no-per-base \\
        --fast-mode \\
        --by 100 \\
        --quantize 0:1:5:20:50: \\
        ${meta.id} \\
        ${cram}
    """

    stub:
    """
    touch ${meta.id}.mosdepth.global.dist.txt ${meta.id}.mosdepth.region.dist.txt ${meta.id}.mosdepth.summary.txt
    for bed in regions quantized; do
        touch ${meta.id}.\$bed.bed.gz ${meta.id}.\$bed.bed.gz.csi
    done
    """
}
