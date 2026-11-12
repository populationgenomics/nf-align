/*
 * CRAM QC: VerifyBamID, samtools stats, Picard metrics and mosdepth coverage of each CRAM, in parallel, one task per
 * tool's image. The CRAMs are newly aligned ones, and with --qc_aligned already-aligned ones too, so each
 * comes with the reference it was written against. An SG is emitted only once all four have succeeded, so its CRAM
 * QC is all or nothing.
 */

include { VERIFYBAMID } from '../../../modules/local/verifybamid'
include { SAMTOOLS_STATS } from '../../../modules/local/samtools_stats'
include { PICARD_METRICS } from '../../../modules/local/picard_metrics'
include { MOSDEPTH } from '../../../modules/local/mosdepth'

workflow CRAM_QC {
    take:
    crams           // [meta, cram, crai, fasta, fai, dict], each CRAM with the reference it was written against
    contamination   // value: [ud, bed, mu], VerifyBamID's resources for the sequencing type
    intervals       // value: genome coverage, or exome evaluation, BED format
    sequencing_type // string: genome | exome

    main:
    def verify = VERIFYBAMID(crams, contamination, sequencing_type)
    def stats = SAMTOOLS_STATS(crams)
    def picard = PICARD_METRICS(crams, intervals, sequencing_type)
    def coverage = MOSDEPTH(crams)

    // keyed by SG ID; join drops an SG missing from any of the four
    def complete = verify.map { meta, selfsm -> [meta.id, meta, [selfsm]] }
        .join(stats.map { meta, samtools -> [meta.id, [samtools]] })
        .join(picard.map { meta, f1, f2, f3, f4, f5, f6 -> [meta.id, [f1, f2, f3, f4, f5, f6]] })
        .join(coverage.map { meta, f1, f2, f3, f4, _f5, f6, _f7 -> [meta.id, [f1, f2, f3, f4, f6]] })
        .map { _id, meta, selfsm, samtools, metrics, depth -> [meta, selfsm + samtools + metrics + depth] }

    emit:
    qc = complete // [meta, [the SG's QC files]], for SGs whose every QC task succeeded
}
