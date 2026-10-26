/*
 * FASTQ / FASTQ.ora / BAM / CRAM -> deduplicated, coordinate-sorted CRAM, with bwa-mem3 and dupblaster.
 * A single unit (one FASTQ pair, or an unsharded CRAM/BAM) aligns, dedups and writes the CRAM in one task.
 * Multiple units (FASTQ lanes, or bazam shards) align separately and are merged before dedup.
 * SGs that were already aligned pass straight through to `all_crams` without adding any tasks.
 */

include { SORT_INDEX } from '../../../modules/local/sort_index'
include { ALIGN_TO_CRAM } from '../../../modules/local/align_to_cram'
include { ALIGN_SHARD } from '../../../modules/local/align_shard'
include { MERGE_TO_CRAM } from '../../../modules/local/merge_to_cram'

workflow ALIGN {
    take:
    to_align        // [meta, sample sheet row] with meta.kind in fastq | fastq_ora | cram | bam
    aligned         // [meta, sample sheet row] with meta.kind == aligned
    reference       // value: [fasta, [fai and bwa-mem3 index files, named <fasta>.<ext>]]
    shards          // int: bazam shards per CRAM/BAM input
    dupblaster_args // string: extra dupblaster arguments, inserted into the command line verbatim

    main:
    def by_kind = to_align.branch { meta, r ->
        unindexed: meta.kind in ['cram', 'bam'] && !r.reads_index
        indexed: meta.kind in ['cram', 'bam']
        fastq: true
    }

    // CRAM/BAM: sort+index if needed, then fan out into bazam shards
    def sorted = SORT_INDEX(by_kind.unindexed.map { meta, r -> [meta, file(r.reads), file(meta.ref), file("${meta.ref}.fai")] })
    def cram_units = by_kind.indexed
        .map { meta, r -> [meta, file(r.reads), file(r.reads_index)] }
        .mix(sorted)
        .combine(channel.of(1..shards))
        .map { meta, reads, index, i ->
            [meta + [units: shards], i, [reads, index, file(meta.ref), file("${meta.ref}.fai")], shards > 1 ? "${i},${shards}" : '']
        }

    // FASTQ / FASTQ.ora: one row per lane, the ORA reference tarball travels with each lane
    def fastq_units = by_kind.fastq.map { meta, r ->
        def reads = [file(r.r1), file(r.r2)] + (meta.kind == 'fastq_ora' ? [file(meta.ref)] : [])
        [meta, r.unit.toInteger(), reads, '']
    }

    def units = fastq_units
        .mix(cram_units)
        .branch { meta, _unit, _reads, _shard ->
            single: meta.units == 1
            sharded: true
        }

    ALIGN_TO_CRAM(units.single, reference, dupblaster_args)

    // groupKey releases each SG as soon as its last shard finishes; sorting keeps -resume hashes stable.
    // An SG that lost a shard to a failed task never completes its group, and groupTuple's default
    // `remainder: false` drops it instead of merging the shards that survived.
    def shard_bams = ALIGN_SHARD(units.sharded, reference)
        .map { meta, bam -> [groupKey(meta, meta.units), bam] }
        .groupTuple(sort: true)
        .map { key, bams -> [key.getGroupTarget(), bams] }
    MERGE_TO_CRAM(shard_bams, reference, dupblaster_args)

    def new_crams = ALIGN_TO_CRAM.out.cram.mix(MERGE_TO_CRAM.out.cram)

    emit:
    cram = new_crams                                                    // [meta, cram, crai], newly aligned only
    markdup = ALIGN_TO_CRAM.out.markdup.mix(MERGE_TO_CRAM.out.markdup) // [meta, metrics]
    all_crams = new_crams.mix(aligned.map { meta, r -> [meta, file(r.reads), file(r.reads_index)] }) // for genotyping
}
