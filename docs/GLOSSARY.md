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
An SG's dataset is the one its sample belongs to, whichever cohort brought it into a run.
_Avoid_: project (except when you specifically mean the Metamist project name, which may carry `-test`)

**Analysis dataset**:
The dataset a run is launched under. It pays for the run, holds its run records, and can see every
dataset its SGs belong to. It need not own any of the SGs.
_Avoid_: launch dataset, home dataset

**Access level**:
Whether a run works in its datasets' production (`full`, `standard`) or test (`test`) namespaces. Test
runs read `<dataset>-test` Metamist projects and write to `-test` buckets. Every SG in a run must be at the
run's access level.
_Avoid_: namespace, environment

**Dataset bucket**:
`gs://cpg-<dataset>-main`, or `gs://cpg-<dataset>-test` at the test access level. Outputs are laid out
under it by fixed convention.

**Cohort**:
A named, curated set of sequencing groups in Metamist, and the unit a run is asked to process. It
belongs to one dataset, but its SGs can belong to others. It can be active or inactive.

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
A Metamist record of an output file (e.g. a `cram`, `gvcf` or `vds`), linked to the SGs or cohorts it
covers. It's trusted when it is active. `completed` is a Metamist status, but it isn't relied on.

**Registration**:
Creating an analysis for a finished output: one that has been published, or that was written directly to
its final path (as a VDS is).
_Avoid_: recording, submission

**ar-guid**:
The analysis-runner GUID that identifies one run. It's stamped on compute, on registered analyses, and
on the run's record directory.
_Avoid_: run ID

### Inputs

**Sample sheet**:
The run's work list, built from Metamist or supplied directly. Each pipeline defines what one row is: in
alignment, an alignment unit and the dataset its SG belongs to; in combining, one SG with its gVCF and the
dataset the SG belongs to.
_Avoid_: manifest, input sheet

**Kind**:
The form of an SG's input in the sample sheet: `fastq`, `fastq_ora`, `cram`, `bam`, or `aligned`. The
pre-flight check can also turn an SG into kind `unregistered`.
_Avoid_: reads_type (the Metamist assay field `kind` is derived from)

**Already-aligned SG**:
An SG with a completed, active CRAM analysis. It's carried through as kind `aligned` (for later
genotyping) but never realigned, and its CRAM is never republished or re-registered. If it has no somalier
fingerprint, one is made (gapfilled), published and registered.

**Unregistered SG**:
An SG with no CRAM analysis, whose CRAM and `.crai` an earlier run published to its destination anyway. It's
carried as kind `unregistered`. It isn't realigned: its CRAM is fingerprinted, and the CRAM and fingerprint
are registered separately (ADR 0011).

**Gapfill**:
Fingerprinting an already-aligned SG that has no completed, active `somalier` analysis.

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
The genome reads are aligned to and output CRAMs are encoded against: the unmasked GATK hg38 FASTA (ADR 0010).
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
aligns alt-aware: bwa-mem3 reads the `.alt` file beside the alignment reference (ADR 0010).

### Outputs and runs

**Published output**:
A file copied from the pipeline into its SG's dataset bucket at its conventional path. A VDS is never
published: it's written directly to its final path.

**Output clash**:
A path that this run would publish, but which already exists. It stops the run before any alignment.

**Run records**:
The per-run directory `alignment_samplesheets/<ar-guid>/` in the analysis dataset bucket, holding the
sample sheet, the index of published CRAMs, and the failed SGs.

**Failed SG**:
An SG that was dropped part-way through a run, at either the `align` or the `register` step. The rest of
the run carries on, and the SG is listed in the run records.

### Combining

**gVCF**:
One SG's single-sample genomic VCF, recorded in Metamist as an analysis of type `gvcf`. Every SG combined
must have exactly one distinct active gVCF. Several active analyses pointing at the same file count as
one.
_Avoid_: VCF (a gVCF also records reference blocks)

**VDS**:
A Hail VariantDataset: many SGs' variant and reference data combined into one dataset. It is built from
gVCFs, and can be built on top of an earlier VDS.
_Avoid_: matrix table, callset

**Base VDS**:
The earlier VDS a combine builds on. SGs in it that are no longer requested are removed, and requested SGs
it lacks are added from their gVCFs.
_Avoid_: previous VDS, prior VDS

**Fresh build**:
A combine that ignores any base VDS and builds the VDS from gVCFs alone.

**Refresh**:
Replacing an SG's data in the VDS with its current gVCF. It's needed when an SG keeps its ID but its gVCF
changes, which the VDS can't detect by itself. A refreshed SG is removed from the base VDS and added back.
_Avoid_: replace (used loosely for a sample getting a new SG)
