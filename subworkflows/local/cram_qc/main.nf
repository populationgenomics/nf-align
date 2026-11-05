/*
 * CRAM QC: VerifyBamID, samtools stats and Picard metrics of each CRAM, in parallel, one task per tool's image. The
 * CRAMs are newly aligned ones, and with --qc_aligned already-aligned ones too (ADR 0015), so each comes with the
 * reference it was written against. An SG is emitted only once all three have succeeded, so its CRAM QC is all or
 * nothing (ADR 0014).
 */

include { VERIFYBAMID } from '../../../modules/local/verifybamid'
include { SAMTOOLS_STATS } from '../../../modules/local/samtools_stats'
include { PICARD_METRICS } from '../../../modules/local/picard_metrics'

workflow CRAM_QC {
    take:
    crams           // [meta, cram, crai, fasta, fai, dict], each CRAM with the reference it was written against
    contamination   // value: [ud, bed, mu], VerifyBamID's resources for the sequencing type
    intervals       // value: genome coverage, or exome evaluation, interval list
    sequencing_type // string: genome | exome

    main:
    def verify = VERIFYBAMID(crams, contamination, sequencing_type)
    def stats = SAMTOOLS_STATS(crams)
    def picard = PICARD_METRICS(crams, intervals, sequencing_type)

    // keyed by SG ID; join drops an SG missing from any of the three
    def complete = verify.map { meta, selfsm -> [meta.id, meta, [selfsm]] }
        .join(stats.map { meta, samtools -> [meta.id, [samtools]] })
        .join(picard.map { meta, f1, f2, f3, f4, f5, f6 -> [meta.id, [f1, f2, f3, f4, f5, f6]] })
        .map { _id, meta, selfsm, samtools, metrics -> [meta, selfsm + samtools + metrics] }

    emit:
    qc = complete // [meta, [the SG's QC files]], for SGs whose every QC task succeeded
}
