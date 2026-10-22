// single unit (one FASTQ pair, or an unsharded CRAM/BAM): align, dedup and write the CRAM in one task
process ALIGN_TO_CRAM {
    tag "${meta.id}"

    input:
    tuple val(meta), val(unit), path(reads, stageAs: 'in/*'), val(shard)
    tuple val(index_name), path(index_files, stageAs: 'bwa_index/*')
    tuple path(fasta), path(fai)
    val dupblaster_args

    output:
    tuple val(meta), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), emit: cram
    tuple val(meta), path("${meta.id}.markduplicates-metrics"), emit: markdup

    script:
    """
    eval "\$(threads.sh ${task.cpus})"
    align_stream.sh ${meta.kind} '${shard}' bwa_index/${index_name} ${meta.id} \$BWA_THREADS \$READ_THREADS ${reads} \\
      | dedup_to_cram.sh ${meta.id} ${fasta} \$SORT_THREADS \$WRITE_THREADS ${dupblaster_args}
    """

    stub:
    """
    touch ${meta.id}.cram ${meta.id}.cram.crai ${meta.id}.markduplicates-metrics
    """
}
