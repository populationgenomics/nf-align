/*
 * Picard metrics of a CRAM, as cpg-flow's CramQcPicardMultiMetrics and CramQcPicardCollectMetrics, in one task
 * so the CRAM is staged once: CollectMultipleMetrics runs alongside CollectWgsMetrics (genome) or
 * CollectHsMetrics (exome), each JVM with 40% of the task's memory.
 *
 * htsjdk's CRAM decoding is single-threaded, and was most of the task's time, so samtools decodes the CRAM
 * instead, multi-threaded, and streams each JVM uncompressed BAM on stdin.
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
    // decode threads per stream, on top of samtools' own
    def threads = Math.max(1, task.cpus.intdiv(2) - 1)
    // a stream has no index, so picard can't do the indexed lookup of the coverage regions it would on the CRAM:
    // samtools does it instead (-M -L: each read overlapping the regions, once). HsMetrics needs every read
    def regions = sequencing_type == 'genome' ? '-M -L intervals.bed' : ''
    // the interval lists' headers carry the md5s of the reference they were made against (the DRAGEN-masked
    // one, for the genome coverage regions), and picard fails on a mismatch. A round trip through BED drops the
    // header, and re-adds it from this CRAM's .dict: the same regions, against the right reference
    def coverage = sequencing_type == 'genome'
        ? """
        decode ${regions} | picard ${heap} CollectWgsMetrics \\
            -I /dev/stdin \\
            -R ${fasta} \\
            -O ${meta.id}.picard-wgs-metrics \\
            --VALIDATION_STRINGENCY SILENT \\
            --INTERVALS intervals.interval_list \\
            --USE_FAST_ALGORITHM true \\
            --READ_LENGTH 250
        """
        : """
        decode | picard ${heap} CollectHsMetrics \\
            -I /dev/stdin \\
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
    # the CRAM, as uncompressed BAM. A reader may close the stream before its end (CollectWgsMetrics stops after
    # its last interval), and samtools' SIGPIPE (141) then isn't a failure: each pipeline's status is picard's
    decode() {
        samtools view -@${threads} -u -T ${fasta} "\$@" ${cram} || [ \$? -eq 141 ]
    }

    picard IntervalListToBed -I ${intervals} -O intervals.bed
    picard BedToIntervalList -I intervals.bed -O intervals.interval_list -SD ${dict}

    decode | picard ${heap} CollectMultipleMetrics \\
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
        -LEVEL SAMPLE &
    multiple=\$!

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
