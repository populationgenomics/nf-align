// one of several units (a FASTQ lane, or a bazam shard): align and name-sort for MERGE_TO_CRAM
process ALIGN_SHARD {
    tag "${meta.id}:${unit}"

    input:
    tuple val(meta), val(unit), path(reads, stageAs: 'in/*'), val(shard)
    tuple val(index_name), path(index_files, stageAs: 'bwa_index/*')

    output:
    tuple val(meta), path("${meta.id}.${unit}.bam")

    script:
    """
    eval "\$(threads.sh ${task.cpus})"
    mkdir -p tmp
    align_stream.sh ${meta.kind} '${shard}' bwa_index/${index_name} ${meta.id} \$BWA_THREADS \$READ_THREADS ${reads} \\
      | samtools sort -n -@\$SORT_THREADS -T tmp/sort -o ${meta.id}.${unit}.bam
    check_alignments.sh ${meta.id}.${unit}.bam
    """

    stub:
    """
    touch ${meta.id}.${unit}.bam
    """
}
