/*
 * Picard metrics of a CRAM, as cpg-flow's CramQcPicardMultiMetrics and CramQcPicardCollectMetrics, in one task
 * so the CRAM is staged once: CollectMultipleMetrics runs alongside CollectWgsMetrics (genome) or
 * CollectHsMetrics (exome), each JVM with 40% of the task's memory.
 *
 * htsjdk's CRAM decoding is single-threaded, and was most of the task's time, so samtools decodes the CRAM
 * instead, once and multi-threaded, and tee streams the uncompressed BAM to both JVMs.
 */
process PICARD_METRICS {
    tag "${meta.id}"

    input:
    // the CRAM's own reference, staged apart so it can't collide with the CRAM's name. Picard finds the .dict
    // beside the FASTA, by name
    tuple val(meta), path(cram), path(crai), path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*'), path(dict, stageAs: 'ref/*')
    // genome: the coverage regions; exome: the evaluation regions, as both bait and target
    path intervals, stageAs: 'intervals/*'
    val sequencing_type

    output:
    tuple val(meta),
        path("${meta.id}.alignment_summary_metrics"),
        path("${meta.id}.base_distribution_by_cycle_metrics"),
        path("${meta.id}.insert_size_metrics"),
        path("${meta.id}.quality_by_cycle_metrics"),
        path("${meta.id}.quality_yield_metrics"),
        path("${meta.id}.picard-${sequencing_type == 'genome' ? 'wgs' : 'hs'}-metrics")

    script:
    def heap = "-Xmx${(task.memory.toMega() * 0.4) as int}m"
    // decode threads, on top of samtools' own. The JVMs take a core each, but a decode thread with nowhere to
    // write blocks on the pipe, so these only run when picard keeps up
    def threads = Math.max(1, task.cpus - 1)
    // the interval lists' headers carry the md5s of the reference they were made against (the DRAGEN-masked
    // one, for the genome coverage regions), and picard fails on a mismatch. A round trip through BED drops the
    // header, and re-adds it from this CRAM's .dict: the same regions, against the right reference.
    // A stream has no index, so picard reads every record and skips those outside the regions itself
    def coverage = sequencing_type == 'genome'
        ? """
        picard ${heap} CollectWgsMetrics \\
            -I /dev/stdin \\
            -R ${fasta} \\
            -O ${meta.id}.picard-wgs-metrics \\
            --VALIDATION_STRINGENCY SILENT \\
            --INTERVALS ${intervals} \\
            --USE_FAST_ALGORITHM true \\
            --READ_LENGTH 250 \\
            < coverage.fifo &
        """
        : """
        picard ${heap} CollectHsMetrics \\
            -I /dev/stdin \\
            -R ${fasta} \\
            -O ${meta.id}.picard-hs-metrics \\
            --VALIDATION_STRINGENCY SILENT \\
            -TI ${intervals} \\
            -BI ${intervals} \\
            -LEVEL null \\
            -LEVEL SAMPLE \\
            -LEVEL LIBRARY \\
            < coverage.fifo &
        """
    """

    # each JVM reads the CRAM, decoded once, as uncompressed BAM from its own FIFO
    mkfifo multiple.fifo coverage.fifo

    picard ${heap} CollectMultipleMetrics \\
        -I /dev/stdin \\
        -R ${fasta} \\
        -O ${meta.id} \\
        -AS True \\
        --VALIDATION_STRINGENCY SILENT \\
        -PROGRAM null \\
        -PROGRAM CollectAlignmentSummaryMetrics \\
        -PROGRAM CollectInsertSizeMetrics \\
        -PROGRAM MeanQualityByCycle \\
        -PROGRAM CollectBaseDistributionByCycle \\
        -PROGRAM CollectQualityYieldMetrics \\
        -LEVEL null \\
        -LEVEL SAMPLE \\
        < multiple.fifo &
    multiple=\$!

    ${coverage.trim()}
    coverage=\$!

    # CollectWgsMetrics closes its FIFO after its last interval: -p (GNU tee) drops a reader that's gone, and
    # carries on feeding the other, rather than dying of SIGPIPE. A reader's failure is caught by its wait
    samtools view -@${threads} -u -T ${fasta} ${cram} | tee -p coverage.fifo > multiple.fifo

    # under set -e, a failed JVM fails the task here
    wait \$multiple
    wait \$coverage
    """

    stub:
    def coverage = sequencing_type == 'genome' ? 'wgs' : 'hs'
    """
    for ext in alignment_summary_metrics base_distribution_by_cycle_metrics insert_size_metrics \\
        quality_by_cycle_metrics quality_yield_metrics picard-${coverage}-metrics; do
        touch ${meta.id}.\$ext
    done
    """
}
