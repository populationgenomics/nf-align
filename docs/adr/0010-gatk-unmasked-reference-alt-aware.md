# Align to the unmasked GATK hg38 FASTA, alt-aware

bwa-mem3 aligns to the GATK `Homo_sapiens_assembly38.fasta`
(`gs://cpg-common-main/rd_references/alignment/unmasked/`), not the DRAGEN-masked FASTA. Its bwa-mem3
index is built once, outside the pipeline, and sits beside it. The alignment is alt-aware: `ref_suffixes`
in `main.nf` includes `alt`, so `<ref_fasta>.alt` is staged beside the index, and bwa-mem3 reads it.
The alt-annotated `.dict` (`ref_fasta_dict`) is staged too, so the `@SQ` lines keep their `AH` tags.

This supersedes [ADR 0006](0006-dragen-masked-reference-alt-unaware.md).

## Consequences

- **Two references across the cohort:** new CRAMs are encoded against a different FASTA from the existing
  CPG CRAMs, which used the masked one. Decoding a CRAM, or analysing new and old CRAMs together, has to use
  each CRAM's own reference. Mixing them in one joint analysis needs a deliberate decision. The old CRAMs
  were also aligned alt-unaware, so they differ in how reads near alt haplotypes are placed, too.
- **The `.alt` and `.dict` travel together:** bwa-mem3 copies the `@SQ` lines of `<baseprefix>.dict`
  (`Homo_sapiens_assembly38.dict`) into the header verbatim. The staged one was made with
  `samtools dict --alt`, so its alt contigs carry `AH:*`. A `.dict` without those tags would make the alt
  contigs lose `AH` in every CRAM header. Swapping either file, or dropping `alt` from `ref_suffixes`,
  changes how reads are placed and needs a new ADR.
- **The stub tests don't cover this:** their placeholder reference has no `.alt` or `.dict`, and a stub
  never reads them. Only real CRAMs show whether alt awareness took effect (RISKS #22).
