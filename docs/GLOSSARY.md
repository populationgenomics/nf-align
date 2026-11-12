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
The product of one sequencing event for a sample, as recorded in Metamist. Its metadata records the reads' type and
locations. An SG is made up of one or more assays combined.

**Analysis**:
A Metamist record of an output file (e.g. a `cram`, `gvcf` or `vds`), linked to the SGs or cohorts it
covers. It's trusted when it is active. `completed` is a Metamist status, but it isn't relied on.

**Registration**:
Creating an analysis for a finished output: one that has been published, or that was written directly to
its final path (as a VDS is).
_Avoid_: recording, submission

**ar-guid**:
The analysis-runner GUID that identifies one run. It's stamped on compute, on registered analyses, and
on the run's record directory. It is also used to aggregate billing records for a NextFlow run.
_Avoid_: run ID

### Inputs

**Sample sheet**:
The run's work list, built from Metamist or supplied directly. Each pipeline defines what one row is: in
alignment, an alignment unit and the dataset its SG belongs to; in genotyping, one SG with its CRAM and
the dataset the SG belongs to; in combining, one SG with its gVCF and the dataset the SG belongs to.
_Avoid_: manifest, input sheet

**Kind**:
The form of an SG's input in the sample sheet: `fastq`, `fastq_ora`, `cram`, `bam`, or `aligned`.
_Avoid_: reads_type (the Metamist assay field `kind` is derived from)

**Already-aligned SG**:
An SG with a completed, active CRAM analysis. It's carried through as kind `aligned` (for later
genotyping) but never realigned, refingerprinted, republished or re-registered, whether or not its CRAM
analysis has a somalier fingerprint (ADR 0012). With `--qc_aligned`, one without CRAM QC is a QC-only SG.

**QC-only SG**:
An SG that gets QC and nothing else. In alignment, it's an already-aligned SG whose CRAM has no `cram_qc`
analysis, in a run with `--qc_aligned true`, and it gets CRAM QC against its CRAM's own reference (ADR 0015). In
genotyping, it's a genotyped SG whose gVCF has no `gvcf_qc` analysis, in a run with `--qc_genotyped true`, and it
gets gVCF QC against its gVCF's own reference.
_Avoid_: QC gapfill SG, re-QC

**Unregistered CRAM**:
A CRAM at an SG's destination with no CRAM analysis, for example because the run that published it failed to
fingerprint or register it. It's an output clash: nf-align only registers files it wrote (ADR 0013).

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
The per-library summary of duplicate marking. It's written by alignment, but it's part of the SG's
CRAM QC.

**Alt contigs**:
The alternative haplotype sequences in hg38. An alt-aware aligner treats them specially. This pipeline
aligns alt-aware: bwa-mem3 reads the `.alt` file beside the alignment reference (ADR 0010).

### Outputs and runs

**Published output**:
A file copied from the pipeline into its SG's dataset bucket at its conventional path. A VDS is never
published: it's written directly to its final path.

**Output clash**:
A path that this run would publish, but which already exists. It stops the run before any alignment, at
`PREFLIGHT`. A resumed run whose pre-flight check passed doesn't check again, so its own earlier attempt's files
are never clashes (ADR 0016).

**Run records**:
The per-run directory in the analysis dataset bucket that holds the sample sheet:
`alignment_samplesheets/<ar-guid>/` for alignment, `genotyping_samplesheets/<ar-guid>/` for genotyping. Neither the CRAMs the run produced (ADR 0019) nor its failed SGs (ADR 0017) are recorded there.

**Failed SG**:
An SG that didn't get through every step of a run: alignment, fingerprinting, registration, CRAM QC or its
registration. The rest of the run carries on, and the SG isn't recorded anywhere but in Seqera's task list;
`-resume` picks it up (ADR 0017). An SG that fails only at QC or its registration still has its CRAM
registered.

**CRAM QC**:
Quality-control metrics computed from one CRAM (newly aligned, or a QC-only SG's), kept together in one folder per SG and
recorded as one `cram_qc` analysis. Its primary file is the CRAM, and every metrics file is a secondary
file. It's all or nothing: an SG either gets every file registered or no `cram_qc` analysis (a QC-only
SG's may lack the duplicate metrics, which aren't this run's). It never
holds up the CRAM's own registration, and it's only registered after the CRAM is.
_Avoid_: align QC, alignment QC, `qc` analysis (a different, dataset-level type)

**CRAM companion**:
A file that belongs to the CRAM itself and is registered on its `cram` analysis: the index and the somalier
fingerprint. QC metrics are not companions.

### Genotyping

**Genotyping**:
Turning one SG's CRAM into its gVCF. It always follows the CRAM's alignment and CRAM QC. There is no
multi-sample genotyping step: combining gVCFs into a VDS is a separate process.
_Avoid_: variant calling, calling (as the name of the step)

**gVCF**:
One SG's variant calls, with the stretches between them summarised as blocks of reference confidence, made
from one CRAM. An SG has at most one active gVCF analysis at a time.
_Avoid_: VCF (a gVCF is not a plain VCF), g.vcf

**Source CRAM**:
The CRAM a gVCF was genotyped from. The gVCF is genotyped against the source CRAM's own reference, and both
are recorded on the gVCF's analysis.

**Genotyped SG**:
An SG that already has an active gVCF analysis. It's never genotyped again, whichever CRAM that gVCF came from.
_Avoid_: called SG

**Genotyping interval**:
One of the sets of genomic regions that an SG's genotyping is split into, each genotyped separately and
then merged back into one gVCF. Every SG in a run is split the same way.
_Avoid_: shard (a read subset, in alignment), chunk, calling interval

**gVCF QC**:
Quality-control metrics computed from one gVCF, kept together in one folder per SG and recorded as one
`gvcf_qc` analysis. Its primary file is the gVCF, and every metrics file is a secondary file. It's all or
nothing, and it's only registered after the gVCF is (a back-filled one, from before genotyping wrote
`bcftools stats`, lacks that file, so it's optional). Aggregating it across SGs is not part of genotyping.
_Avoid_: variant calling QC, VCF QC
