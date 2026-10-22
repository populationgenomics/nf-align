# nf-align

Alignment of short-read sequencing data held in Metamist. Each sequencing group's reads are turned into
one deduplicated, coordinate-sorted CRAM, which is then recorded back in Metamist.

## Language

### Metamist and CPG

**Metamist**:
CPG's sample-metadata service, and the system of record for projects, cohorts, sequencing groups, assays
and analyses.

**Dataset**:
A Metamist project, named without any `-test` suffix. It owns the buckets that outputs are written to.
_Avoid_: project (except when you specifically mean the Metamist project name, which may carry `-test`)

**Access level**:
Whether a run works in a dataset's production (`full`, `standard`) or test (`test`) namespace. Test runs
read the `<dataset>-test` Metamist project and write to the `-test` bucket.
_Avoid_: namespace, environment

**Dataset bucket**:
`gs://cpg-<dataset>-main`, or `gs://cpg-<dataset>-test` at the test access level. Outputs are laid out
under it by fixed convention.

**Cohort**:
A named, curated set of sequencing groups in Metamist, and the unit a run is asked to process. It can be
active or inactive.

**Sequencing group (SG)**:
All of the sequencing data for one sample, of one sequencing type and technology. It is the unit of
alignment: one SG gives one CRAM. Its ID (e.g. `CPG123456`) names its outputs and its read group.
_Avoid_: sample (a sample can have several SGs)

**Sequencing type**:
`genome` (WGS) or `exome` (WES). Every run handles exactly one.
_Avoid_: WGS/WES in parameters and paths

**Technology**:
How the data was sequenced. Only `short-read` is in scope.

**Assay**:
One sequencing event for an SG, as recorded in Metamist. Its metadata records the reads' type and
locations. An SG can have several assays, and their reads are combined.

**Analysis**:
A Metamist record of an output file (here, a `cram`), linked to the SGs it covers. It's only trusted
when it is completed and active.

**Registration**:
Creating an analysis for an output that has been published.
_Avoid_: recording, submission

**ar-guid**:
The analysis-runner GUID that identifies one run. It's stamped on compute, on registered analyses, and
on the run's record directory.
_Avoid_: run ID

### Inputs

**Sample sheet**:
The run's work list: one row per alignment unit. It's built from cohorts or SG IDs, or supplied directly.
_Avoid_: manifest, input sheet

**Kind**:
The form of an SG's input in the sample sheet: `fastq`, `fastq_ora`, `cram`, `bam`, or `aligned`.
_Avoid_: reads_type (the Metamist assay field `kind` is derived from)

**Already-aligned SG**:
An SG with a completed, active CRAM analysis. It's carried through as kind `aligned` (for later
genotyping) but never realigned, republished or re-registered.

**Lane**:
One FASTQ R1/R2 pair, contributed by one assay.
_Avoid_: read pair (a single template), unit

**FASTQ.ora**:
FASTQ compressed with Illumina's ORA format. It needs an ORA reference to decompress.

**Input reference**:
The reference an input CRAM was encoded against. It's needed to decode that CRAM, and it can differ from
the alignment reference.
_Avoid_: old reference

### Alignment

**Alignment reference**:
The genome reads are aligned to and output CRAMs are encoded against: the DRAGEN-masked hg38 FASTA.
_Avoid_: reference (unqualified, when the input reference could be meant)

**Alignment unit**:
The smallest piece of an SG that is aligned on its own: one lane, one shard, or a whole unsharded
CRAM/BAM.
_Avoid_: chunk, part

**Shard**:
A subset of a CRAM/BAM's read pairs, chosen by a hash of the read name, so both mates always land in the
same shard. Shards don't correspond to genomic regions.
_Avoid_: region, interval, split

**Read group**:
The tag that ties each read to its source. Here there's one per SG, with ID and SM both set to the SG ID.

**Query-grouped**:
Read order in which all records with the same read name sit next to each other. The relative order of
different names doesn't matter. Aligner output and name-sorted merges are both query-grouped, and
duplicate marking requires it.
_Avoid_: RG-ordered, name-sorted (a stricter property)

**Duplicate marking**:
Flagging read pairs that are copies of the same original DNA fragment, so later analyses can ignore
them. It's done once per SG, across all of its units.
_Avoid_: dedup, deduplication (reads are flagged, not removed)

**Sequencing duplicate**:
A duplicate created by the sequencer (optical or pad-hopping) rather than by library amplification. It's
told apart using coordinates parsed from Illumina read names.
_Avoid_: optical duplicate (a subset)

**Duplicate metrics**:
The per-library summary of duplicate marking, kept alongside each CRAM as
`<sg>.markduplicates-metrics`.

**Alt contigs**:
The alternative haplotype sequences in hg38. An alt-aware aligner treats them specially. This pipeline
aligns without alt awareness.

### Outputs and runs

**Published output**:
A file copied from the pipeline into the dataset bucket at its conventional path. Only published
outputs are registered.

**Output clash**:
A path that this run would publish, but which already exists. It stops the run before any alignment.

**Run records**:
The per-run directory `alignment_samplesheets/<ar-guid>/`, holding the sample sheet, the index of
published CRAMs, and the failed SGs.

**Failed SG**:
An SG that was dropped part-way through a run, at either the `align` or the `register` step. The rest of
the run carries on, and the SG is listed in the run records.
