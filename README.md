# nf-align

Alignment repository built in NextFlow.

Takes Metamist cohorts (or a list of sequencing group IDs), which can span several datasets, and produces a
deduplicated, coordinate-sorted CRAM 3.0 for each short-read sequencing group, whether its input is FASTQ,
FASTQ.ora, BAM or CRAM. It aligns with bwa-mem3 and marks duplicates with dupblaster. Each SG's CRAM is
published into its own dataset's bucket (`cram/<SG>.cram`, or `exome/cram/<SG>.cram`) and registered in that
dataset's Metamist project. The run records go to the analysis dataset (`--dataset`). It runs on Google Batch.

```bash
nextflow run main.nf \
    --dataset <analysis dataset> --sequencing_type genome --arGuid <guid> --cohorts COH1,COH2 \
    --align_container <image> --cpg_container <image>

# offline wiring check (outputs under results/<dataset>/), and the tests
nextflow run main.nf -profile test -stub --dataset ds --sequencing_type genome \
    --samplesheet sheet.tsv --ref_fasta refs/ref.fasta
uv sync && uv run pytest
```

Documentation:
- [docs/DESIGN.md](docs/DESIGN.md): design brief, flow, outputs, resources, status and next steps
- [docs/GLOSSARY.md](docs/GLOSSARY.md): domain terms
- [docs/adr/](docs/adr/): decision records
- [docs/RISKS.md](docs/RISKS.md): what may bite later
- [docs/AMBIGUOUS-CHOICES.md](docs/AMBIGUOUS-CHOICES.md): judgement calls, including ones made during implementation that still need review
