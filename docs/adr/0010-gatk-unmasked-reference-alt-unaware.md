# Align to the unmasked GATK hg38 FASTA, without alt awareness

bwa-mem3 aligns to the GATK `Homo_sapiens_assembly38.fasta`
(`gs://cpg-common-main/rd_references/alignment/unmasked/`), not the DRAGEN-masked FASTA. Its bwa-mem3
index is built once, outside the pipeline, and sits beside it. As before, no `.alt` file is staged
(`ref_suffixes` in `main.nf` doesn't include `alt`), so the alignment is alt-unaware.

This supersedes [ADR 0006](0006-dragen-masked-reference-alt-unaware.md).

## Consequences

- **Two references across the cohort:** new CRAMs are encoded against a different FASTA from the existing
  CPG CRAMs, which used the masked one. Decoding a CRAM, or analysing new and old CRAMs together, has to use
  each CRAM's own reference. Mixing them in one joint analysis needs a deliberate decision.
- **Alt contigs without alt awareness:** the FASTA includes the alt contigs, and bwa-mem3 treats them as
  ordinary sequence. Reads in regions that have an alt haplotype can map equally well to both, and then get
  MAPQ 0. Check this in the first real CRAMs (RISKS #22).
- **The `.alt` must stay unstaged:** if `alt` is ever added to `ref_suffixes`, bwa-mem3 switches to
  alt-aware mapping. That also needs the `.dict` regenerated with `samtools dict --alt`, as described in
  ADR 0006's considered options.
