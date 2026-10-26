// dedup runs after the merge, as dupblaster must see every shard to catch cross-shard duplicates
process MERGE_TO_CRAM {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bams, stageAs: 'shards/*')
    tuple path(fasta), path(ref_files)
    val dupblaster_args

    output:
    tuple val(meta), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), emit: cram
    tuple val(meta), path("${meta.id}.markduplicates-metrics"), emit: markdup

    script:
    """
    eval "\$(threads.sh ${task.cpus} --no-aligner)"
    samtools merge -n -u -@\$READ_THREADS - ${bams} \\
      | dedup_to_cram.sh ${meta.id} ${fasta} \$SORT_THREADS \$WRITE_THREADS ${dupblaster_args}
    """

    stub:
    """
    touch ${meta.id}.cram ${meta.id}.cram.crai ${meta.id}.markduplicates-metrics
    """
}
