/*
 * Picard metrics of a CRAM, as cpg-flow's CramQcPicardMultiMetrics and CramQcPicardCollectMetrics, in one task
 * so the CRAM is staged once: CollectMultipleMetrics runs alongside CollectWgsMetrics (genome) or
 * CollectHsMetrics (exome), each JVM with 40% of the task's memory.
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
    // the interval lists' headers carry the md5s of the reference they were made against (the DRAGEN-masked
    // one, for the genome coverage regions), and picard fails on a mismatch. A round trip through BED drops the
    // header, and re-adds it from this CRAM's .dict: the same regions, against the right reference
    def coverage = sequencing_type == 'genome'
        ? """
        picard ${heap} CollectWgsMetrics \\
            -I ${cram} \\
            -R ${fasta} \\
            -O ${meta.id}.picard-wgs-metrics \\
            --VALIDATION_STRINGENCY SILENT \\
            --INTERVALS intervals.interval_list \\
            --USE_FAST_ALGORITHM true \\
            --READ_LENGTH 250
        """
        : """
        picard ${heap} CollectHsMetrics \\
            -I ${cram} \\
            -R ${fasta} \\
            -O ${meta.id}.picard-hs-metrics \\
            --VALIDATION_STRINGENCY SILENT \\
            -TI intervals.interval_list \\
            -BI intervals.interval_list \\
            -LEVEL null \\
            -LEVEL SAMPLE \\
            -LEVEL LIBRARY
        """
    """
    picard ${heap} CollectMultipleMetrics \\
        -I ${cram} \\
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
        -LEVEL SAMPLE &
    multiple=\$!

    picard IntervalListToBed -I ${intervals} -O intervals.bed
    picard BedToIntervalList -I intervals.bed -O intervals.interval_list -SD ${dict}
    ${coverage}

    # under set -e, a failed CollectMultipleMetrics fails the task here
    wait \$multiple
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
