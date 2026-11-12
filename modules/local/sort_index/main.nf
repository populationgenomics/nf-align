// an unindexed CRAM/BAM can't be sharded by bazam, so coordinate-sort and index it first
// we anticipate that most BAM/CRAM files as alignment input will arrive with indexes, but safeguard in case
process SORT_INDEX {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads, stageAs: 'in/*'), path(ref, stageAs: 'in/*'), path(fai, stageAs: 'in/*')

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai")

    script:
    """
    mkdir -p tmp
    samtools sort --reference ${ref} -@${task.cpus - 1} -T tmp/sort -o ${meta.id}.sorted.bam ${reads}
    samtools index -@${task.cpus - 1} ${meta.id}.sorted.bam
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    """
}
