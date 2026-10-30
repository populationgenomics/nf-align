---
status: superseded by ADR-0010
---

# Align to the DRAGEN-masked hg38 FASTA, without alt awareness

bwa-mem3 aligns to `Homo_sapiens_assembly38_masked.fasta`, the reference the existing CPG CRAMs were made
with. Its index is built once, outside the pipeline, and no `.alt` file is supplied. Keeping the same
reference means new and existing CRAMs decode against one FASTA and sit side by side in joint analyses.
Changing this later means realigning, or at least re-encoding, every CRAM.

## Considered Options

- **The GATK `Homo_sapiens_assembly38.fasta`, aligned alt-aware** (its `.64.alt` copied to
  `<prefix>.alt`). This would also need the `.dict` regenerated with `samtools dict --alt`, otherwise
  bwa-mem3 drops the `AH` tags from `@SQ`. It was rejected so that the existing cohort stays on one
  reference.
