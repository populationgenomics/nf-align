# nf-align: design brief

Aligns the short-read sequencing groups of a set of Metamist cohorts, or a list of sequencing groups,
which can come from several datasets. The output for each is a deduplicated, coordinate-sorted CRAM with a
somalier fingerprint and a set of CRAM QC metrics, published to and registered in the SG's own dataset. On
request, SGs that are already aligned get CRAM QC alone. It's a Nextflow (≥ 26.04, strict syntax) pipeline
that runs on Google Batch against data in GCS, and it can be tested offline with stub runs.

The pipeline adapts the alignment half of the Nextflow port in `cpg-flow-align-genotype/nextflow/`.
It aligns with **bwa-mem3** (that port used DRAGMAP), and adds an in-pipeline sample sheet step,
fingerprinting, CRAM QC, Metamist registration and a pre-flight output check. Genotyping is left out for
now, but it is the next addition, and the shapes below are chosen so it can be wired in from `main.nf` alone.

Terms are defined in [GLOSSARY.md](GLOSSARY.md). [RISKS.md](RISKS.md) lists what may bite later, and
[AMBIGUOUS-CHOICES.md](AMBIGUOUS-CHOICES.md) lists the calls that could reasonably have gone the other way. This document
was built up in phases through creation of individual Architecture Design Records (ADRs), but during
development a number of these choices were scrapped, replaced, or superseded in some way. This design doc
represents the summary of those documents, but this codebase does not retain the individual files.

## Scope of a run

One run is launched under one **analysis dataset** and covers **one sequencing type**. Its SGs can come
from any number of datasets. Each SG is published to its own dataset's bucket and registered in its own
dataset's Metamist project.

An SG's dataset is its sample's Metamist project, with any `-test` removed. The cohort it came through
doesn't matter: a cohort in `alpha` that holds an SG from `beta` sends that SG's CRAM to `beta`'s bucket.

| param | meaning |
|---|---|
| `dataset` | the analysis dataset (no `-test` suffix): the run is launched under it, it pays for the run, and it holds the run records. Its service account must be able to see every SG's dataset and write to their buckets and Metamist projects |
| `sequencing_type` | `genome` or `exome` |
| `accessLevel` | `full` / `standard` (equivalent), or `test`. Applies to the whole run: every SG's Metamist project must be `<ds>` (full) or `<ds>-test` (test), and its bucket is `gs://cpg-<ds>-main` or `-test` |
| `arGuid` | the analysis-runner GUID, used in resource labels, the per-run record directory and the registered analyses. Mandatory |
| `cohorts` **or** `sgids` **or** `samplesheet` | exactly one: Metamist cohort IDs, Metamist SG IDs (queried directly), or a prepared sheet that skips the Metamist step |
| `output_root` | null by default, so each dataset's outputs go to its bucket. If set (the `test` profile sets `results`), they go to `<output_root>/<ds>/` instead, with the same layout below that |
| `ref_fasta` | the alignment reference. Its `.fai`, bwa-mem3 index and `.alt` (`<ref_fasta>.{fai,bwt.2bit.64,amb,ann,pac,alt}`) must sit beside it, and `ref_fasta_dict` is its alt-annotated `.dict`. Aligned alt-aware |
| `aligned_ref_fasta` | the reference of an already-aligned CRAM whose analysis doesn't record one: the DRAGEN-masked FASTA the CRAMs made before nf-align used. `BUILD_SAMPLESHEET` writes it into such a CRAM's row. Its `.fai` must sit beside it, and, for QC-only SGs, its `<name>.dict` |
| `align_container`, `cpg_container`, `somalier_container` | the images, with defaults in `nextflow.config`. The run fails at launch if any is unset, in every profile |

The defaults are in `nextflow.config`, which `conf/base.config` reads. `main.nf`'s `params {}` block declares
each one's type, so Nextflow converts command-line values (`--cram_qc false` is a Boolean, not the string
`'false'`) and fails the run at launch if a param that isn't optional has no value or the wrong type.
`validateParams()` checks only what the types can't: allowed values, and params needed together.

Switches and tuning, all with defaults in `nextflow.config`:

| param                                                                                                                                                                     | default        | meaning                                                                                                                                         |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| `register`                                                                                                                                                                | `true`         | `false` skips all Metamist registration (the `test` profile sets it)                                                                            |
| `cram_qc`                                                                                                                                                                 | `true`         | `false` skips CRAM QC and the `cram_qc` analysis. The duplicate metrics are still published                                                     |
| `qc_aligned`                                                                                                                                                              | `false`        | also QC already-aligned SGs whose CRAM has no `cram_qc` analysis. Needs `cram_qc` |
| `register_wait_minutes`                                                                                                                                                   | 30             | how long registration polls for the published files                                                                                             |
| `realignment_shards`                                                                                                                                                      | 0              | bazam shards per CRAM/BAM input; 0 means 6 for genomes and 1 for exomes                                                                         |
| `dupblaster_extra_args`                                                                                                                                                   | `''`           | appended verbatim to the dupblaster command (see Running)                                                                                       |
| `somalier_sites`                                                                                                                                                          | hg38 sites VCF | where somalier fingerprints each CRAM                                                                                                           |
| `verifybamid_container`, `picard_container`, `mosdepth_container`, `{genome,exome}_contam_{ud,bed,mu}`, `genome_coverage_interval_list`, `exome_evaluation_interval_list` | cpg-flow's     | CRAM QC's images and references. With `cram_qc` on, the run fails at launch if the images or the sequencing type's references are unset         |

There is no `outputDir`: every destination is worked out per SG (see Outputs). The launcher supplies
`-work-dir` and `google.project`; the config never sets them.

**Two references across the cohort.** New CRAMs are written against `ref_fasta` (unmasked GATK hg38,
alt-aware). CRAMs made before nf-align were written against the DRAGEN-masked FASTA, alt-unaware. A CRAM
is only ever decoded with its own reference, which nf-align records as `ref_fasta` in each `cram`
analysis's meta. Analysing new and old CRAMs together in one joint analysis needs a deliberate decision.

## Flow

```
run-records check ─ --cohorts | --sgids ─ BUILD_SAMPLESHEET ─┐   (cached on -resume)
                                                             ▼
sample sheet ─ parse rows, per-SG destinations into meta ─ PREFLIGHT (every row) ─┬─ aligned: nothing more
                                                                                  ├─ QC-only (--qc_aligned)
                                                                                  └─ to align

   to align: ALIGN ─ publishDir: CRAM, .crai, duplicate metrics → SG's dataset bucket
               ├─ SOMALIER_EXTRACT ─ publishDir: .somalier
               │       └─ REGISTER_CRAM (one cram analysis: CRAM + .crai + .somalier) ─┐
               └─ CRAM_QC ─ publishDir: QC metrics ────────────────────────────────────┴─ REGISTER_CRAM_QC (one cram_qc analysis)

   QC-only:  CRAM_QC (the registered CRAM, its own reference) ─ publishDir: QC metrics ─ REGISTER_CRAM_QC
```

### 1. Sample sheet (`BUILD_SAMPLESHEET`, `bin/build_samplesheet.py`)

- One GraphQL query, using the top-level `cohorts(id:)` or `sequencingGroups(id:)` lookup, so it isn't
  limited to one project. It fetches cohort status and the cohorts' SGs (`--cohorts`), or the named SGs
  directly (`--sgids`), with each SG's project, assay meta and completed, active `cram` and `cram_qc`
  analyses. The script takes `--access-level`, not a project.
- Cohorts that are missing (or not visible to the run's credentials), inactive or empty fail the run. So
  does an SG whose project doesn't match the access level: a production project in a test run, or a
  `-test` project in a full run. Every offender is reported in one error.
- `--cohorts`: SGs are de-duplicated across cohorts. SGs that aren't short-read, or aren't
  `sequencing_type`, are skipped with a warning.
- `--sgids`: each SG was asked for by name, so one that's missing (or not visible), or isn't short-read
  `sequencing_type`, fails the run instead of being skipped.
- Read locations are used **exactly as Metamist stores them**.
- Each SG gets one row per alignment unit, with a `kind` (glossary): `fastq` and `fastq_ora` (one row
  per lane), `cram` or `bam` (one row), or `aligned` (it already has a CRAM analysis; one row pointing at
  that CRAM). An `aligned` row's `somalier` is that CRAM analysis's `somalier` secondary file, or empty if it
  has none. It's for reference only: an aligned SG is never fingerprinted.
  Its `cram_qc` is the ID of its CRAM's latest completed, active `cram_qc` analysis (one whose primary file, `outputs.path`, is that
  same CRAM), or empty. With `--qc_aligned`, an empty one makes the SG QC-only.
- Each row records the SG's `dataset`: its project, with any `-test` removed.
- A run with no in-scope SG fails.
- It logs how many SGs are already aligned, and how many of those lack a fingerprint or CRAM QC.
- The sheet is published to `alignment_samplesheets/<arGuid>/samplesheet.tsv` in the analysis dataset
  bucket.
- Every new run builds its sheet: there is no way to hand one in. A resumed run takes it from the cache, so
  it carries on with the SGs its first attempt started with, and doesn't see what that attempt registered.

Columns: `sg_id, dataset, kind, n_units, unit, r1, r2, reads, reads_index, reads_ref, somalier, cram_qc`. For a CRAM,
`reads_ref` is the reference it was written against, and is required: a CRAM is never decoded with `ref_fasta`.
For ORA, it's the ORA reference tarball. For an `aligned` row, it's the `ref_fasta` its analysis's meta records
(nf-align records it), or else `aligned_ref_fasta` (default: the DRAGEN-masked FASTA the CRAMs made before
nf-align used), which `main.nf` passes to the script. So `main.nf` never has to guess: an `aligned`
row without `reads_ref` fails the run like a `cram` row without one, though the builder never writes one. `main.nf` doesn't read the `somalier` column.

### 2. Pre-flight

The pipeline never overwrites a file. CPG's data protection rules make this difficult/impossible, so it is not attmpted.
`publishDir` uses `overwrite: false`, which skips an existing target silently, so clashes are caught
up front:

- **Run records**: before the sample sheet is built, the run fails if
  `alignment_samplesheets/<arGuid>/` in the analysis dataset bucket already holds a `samplesheet.tsv`.
  This guards against an `arGuid` reused by hand while testing. It's skipped
  on `-resume`.
- **Outputs** (`PREFLIGHT`): before anything is aligned, every CRAM, `.crai`, `.somalier` and QC metrics URI
  that the run would publish is checked in its SG's dataset bucket. If any exist, the run fails and lists them all. `aligned` rows
  are not checked: they publish nothing, and an existing CRAM is the expected state. An unregistered CRAM is a
  clash like any other, because the run only registers files it wrote. A QC-only SG's QC tools' files are
  checked the same way, so QC files left by an earlier, failed attempt stop the run until they're deleted.

`PREFLIGHT` is a native (`exec:`) process, not a function, so that Nextflow caches it. Its input is every
sample sheet row and every URI the run would publish, and it holds every row back until the check has passed. A resumed run whose first
attempt passed the check takes it from the cache, so the files that attempt published aren't clashes, and
`publishDir` skips them when it republishes the cached tasks. A check that failed isn't cached, so a resumed
run checks again. If the rows or destinations change between attempts, the check runs again and finds the
earlier attempt's files: launch a fresh run.

Two runs launched close together for the same SG can both pass the check. The second run's publish is
then skipped silently, and it registers the first run's CRAM under its own `ar-guid`. That gap is accepted.

Before this, at launch, `main.nf` validates params and fails on anything missing or malformed (see Scope
of a run). A sheet row with an unknown `kind`, a `cram` or `aligned` row without `reads_ref`, a row without
`dataset`, or
a QC-only SG whose reference isn't a `.fasta`, `.fa` or `.fna` (so Picard can't find its `.dict`) also fails
the run.

### 3. ALIGN subworkflow

| input | path |
|---|---|
| CRAM/BAM without an index | `SORT_INDEX`, then sharded like an indexed input |
| CRAM/BAM | `realignment_shards` bazam shards, defaulting to 6 for genomes and 1 for exomes |
| FASTQ / FASTQ.ora | one unit per lane; ORA files are streamed through `orad` into bwa-mem3 via FIFOs |

- **One unit** goes to `ALIGN_TO_CRAM`: `align_stream.sh | dedup_to_cram.sh` in a single task.
- **Several units** go to `ALIGN_SHARD` (aligned, then name-sorted to BAM). They are then grouped per SG
  with `groupKey`, and `MERGE_TO_CRAM` runs `samtools merge -n | dedup_to_cram.sh`. Duplicates are
  marked once per SG, after the merge.
- **Alignment:** `bwa-mem3 mem -t <⅔ cpus> -K 100000000 --bam -R '@RG\tID:<sg>\tSM:<sg>\tPL:ILLUMINA'`.
  CRAM/BAM input is piped in interleaved from bazam (`-p -`).
- **`dedup_to_cram.sh`:** `dupblaster | samtools sort -u | samtools view -O cram,version=3.0 --write-index`,
  then the dupblaster duplicate metrics are renamed to `<sg>.markduplicates-metrics`.
- **CPU split:** `bin/threads.sh` is the only place the per-tool ratios live:

  | CPUs | bwa-mem3 | bazam / merge input | sort while aligning | final CRAM write |
  |---|---|---|---|---|
  | 16 (align) | 10 | 3 | 2 | 15 |
  | 8 (merge, `--no-aligner`) | – | 4 | 3 | 7 |

- **Guards:** every script runs under `set -euo pipefail`. `check_alignments.sh` runs `samtools quickcheck`
  and requires at least one record, on every shard BAM and every final CRAM.
- **Reference:** `ref_fasta` is staged with its `.fai`, bwa-mem3 index, `.alt` and alt-annotated `.dict`
  side by side, so the FASTA's own name is the index prefix and the `@SQ` lines keep their `AH` tags.
  Swapping the `.alt` or `.dict`, or dropping `alt` from `ref_suffixes` in `main.nf`, changes how reads are
  placed, and needs a deliberate decision recorded here.
- **Takes** only the rows to align, and **emits** `cram`: `[meta, cram, crai]` for each newly aligned SG.
  Genotyping will need every SG's CRAM: mix `ALIGN.out.cram` with the already-aligned rows in `main.nf`, whose
  `meta.cram` is the registered CRAM.

### 4. Outputs

Everything is published with `publishDir`, `mode: 'copy'` and `overwrite: false`. There is no
`output {}` block: it can only publish under a single `outputDir`
([clarification](clarifications/workflow-outputs-single-bucket.md)).

When `main.nf` parses a sample sheet row, it works out the SG's destinations and puts them in `meta`
next to the rest of its inputs:

- `dataset`, and its Metamist project (`<ds>` or `<ds>-test`)
- `cramDir` and `qcDir`: the two folders it publishes to, under `gs://cpg-<ds>-{main|test}` or
  `<output_root>/<ds>`
- `cram`: where its new CRAM is published, or an already-aligned SG's registered CRAM
- `ref`: the reference that CRAM is (or will be) written against. The input reads' own reference stays in the
  sample sheet row, which only `ALIGN` reads

The `publishDir` settings in config only apply these values: the CRAM and its companions go to `cramDir`, and
every QC file to `qcDir`. `outputs(meta)` in `main.nf` lists every file the SG publishes, by its key in the
analysis that registers it, and the pre-flight check and both registrations read that, so publishing and
registration can't drift apart. `main.nf` stays the one place the layout is defined:

| target | path under the SG's output root (`exome/` prefix for exomes) |
|---|---|
| CRAM + index | `[exome/]cram/<sg>.cram`, plus `.crai` |
| somalier fingerprint | `[exome/]cram/<sg>.cram.somalier` (cpg-flow's `CramPath.somalier_path`) |
| duplicate metrics | `[exome/]cram_qc/<sg>/<sg>.markduplicates-metrics`, even with `--cram_qc false` |
| CRAM QC | `[exome/]cram_qc/<sg>/`: `<sg>.verify-bamid.selfSM`, `<sg>.samtools-stats`, `<sg>.{alignment_summary,base_distribution_by_cycle,insert_size,quality_by_cycle,quality_yield}_metrics`, `<sg>.picard-{wgs,hs}-metrics`, `<sg>.mosdepth.{global.dist,region.dist,summary}.txt`, `<sg>.{regions,quantized}.bed.gz`, and their `.csi`, which are published but not registered |

The run records are the sample sheet, `samplesheet.tsv`, published by `BUILD_SAMPLESHEET` to
`alignment_samplesheets/<arGuid>/` under the analysis dataset's output root. There is no index of the CRAMs a
run produced: each registered analysis records its `ar-guid`, and Seqera lists the tasks. There is no record of failed SGs either.

Only newly aligned SGs are published, apart from a QC-only SG's QC tools' files.

Each per-SG task's `dataset` resource label is `meta.dataset`. `BUILD_SAMPLESHEET` and `PREFLIGHT` are
labelled with the analysis dataset. What the label does to billing depends on the GCP project the run is in
(the CPG billing aggregator, `cpg-infrastructure/cpg_infra/billing_aggregator/aggregate/`):

- **Analysis dataset `seqr`**: `seqr.py` charges each labelled row in full to `topic=<label>`, and adds a
  matching credit on `seqr`. Each SG's compute is billed to its own dataset, and the run-level tasks net
  to zero on `seqr`. A view that leaves out the `-credit` rows shows the same cost under `seqr` and under
  the SG's dataset. Without the label, the cost would be spread across all seqr datasets by size.
- **Any other analysis dataset**: `gcp.py` takes `topic` from the GCP project alone, so the whole run is
  billed to the analysis dataset. The label only fills the `dataset` column of the extended billing view,
  for breaking the run's cost down by SG dataset.

### 5. Fingerprint (`SOMALIER_EXTRACT`)

Every newly aligned CRAM is fingerprinted with `somalier extract` at `somalier_sites`, against `ref_fasta`
(only the FASTA and `.fai` are staged). The `.somalier` is published beside the CRAM and registered as a
secondary file of the CRAM's own analysis, so a CRAM whose fingerprint failed isn't registered.

Already-aligned SGs are never fingerprinted (no **gapfill**): their `cram` analysis already exists, and
`cpg_utils.metamist_registration` can only create analyses, not add a secondary file to one. CRAMs
registered before nf-align fingerprinted have no fingerprint, and this pipeline won't make one.

### 6. CRAM QC (`CRAM_QC`)

When `cram_qc` is true (the default), every newly aligned CRAM is QC'd alongside its fingerprint, as
cpg-flow's `CramQc*` stages did, plus mosdepth
coverage. There are four tasks, one per image, each staging the CRAM with the FASTA, `.fai` and `.dict` of the reference it was
written against (`ref_fasta` and `ref_fasta_dict` for a new CRAM):

- `VERIFYBAMID`: VerifyBamID2 against the sequencing type's 1000G contamination resources
  (`<type>_contam_{ud,bed,mu}`), capped at depth 1000 for exomes.
- `SAMTOOLS_STATS`: `samtools stats`, in the align image.
- `PICARD_METRICS`: CollectMultipleMetrics, run alongside CollectWgsMetrics over
  `genome_coverage_interval_list` (genome), or CollectHsMetrics with `exome_evaluation_interval_list` as
  both bait and target (exome). The interval list is passed through BED and back with the `.dict`, so its
  header matches the alignment reference whichever hg38 build it was made against.
- `MOSDEPTH`: `mosdepth --no-per-base --fast-mode --by 100 --quantize 0:1:5:20:50:`, for both sequencing
  types: the depth distribution, mean depth in 100 bp windows, and the genome split into depth bands. Its
  own task, as mosdepth reads the CRAM through its index and can't share Picard's decoded stream.
  `--fast-mode` ignores CIGARs and mate overlaps: deletions count as covered, and overlapping mates twice.

An SG comes out of `CRAM_QC` only if all four succeeded, so CRAM QC is all or nothing. A QC failure never
touches the CRAM's own registration. All of an SG's metrics, including the duplicate metrics from
alignment, go in one folder, `[exome/]cram_qc/<sg>/`, rather than cpg-flow's per-tool `qc/<tool>/` folders.
The file names inside it are cpg-flow's, so MultiQC's filename patterns still match, but path-based
consumers written for cpg-flow (such as its `CramMultiQC` stage) need to read the `cram_qc` analysis or the
new folder.

Already-aligned SGs aren't QC'd unless `--qc_aligned true` is set. Then each one whose sample sheet row has an
empty `cram_qc` is a **QC-only SG**: its registered CRAM goes through `CRAM_QC` against the reference it was
written with (its row's `reads_ref`, with the `<name>.dict` beside the FASTA), and nothing else
runs for it. Only the four tools' files are published, and the pre-flight check covers them. A run whose SGs
are all aligned does only QC.

### 7. Registration (`REGISTER_CRAM`, `REGISTER_CRAM_QC`)

Both are the `REGISTER_ANALYSIS` process, included twice, and neither runs with `--register false`. A run
only registers files it wrote, with one exception: a QC-only SG's duplicate metrics, below. `REGISTER_CRAM`
runs once per newly aligned SG, after its CRAM is fingerprinted, in the SG's own Metamist project. Each CRAM
gets a single `cram` analysis, with its index and fingerprint as secondary files. A CRAM whose fingerprint
failed isn't registered. The URIs come from `meta`, and the task polls for up to `register_wait_minutes` until the CRAM,
`.crai` and `.somalier` all exist, because publishing is asynchronous. It then runs:

```
python3 -m cpg_utils.metamist_registration --project <SG's dataset>[-test] --output <cram> --type cram \
    --secondary index=<cram>.crai somalier=<cram>.somalier --sgs <sg> \
    --meta sequencing_type=<type> ar-guid=<arGuid> ref_fasta=<CRAM's ref>
```

The fingerprint is found at `outputs.secondary_files.somalier.path` on the CRAM analysis.

`REGISTER_CRAM_QC` runs once an SG's CRAM QC is complete *and* its `cram` analysis is registered (by
`REGISTER_CRAM`, or already, for a QC-only SG), so a `cram_qc` analysis never points at an unregistered CRAM. The CRAM is its primary file again, and every QC
file except mosdepth's BED indexes, including the duplicate metrics, is a secondary file:

```
python3 -m cpg_utils.metamist_registration --project <SG's dataset>[-test] --output <cram> --type cram_qc \
    --secondary verifybamid=<dir>/<sg>.verify-bamid.selfSM samtools_stats=<dir>/<sg>.samtools-stats \
        alignment_summary_metrics=... base_distribution_by_cycle_metrics=... insert_size_metrics=... \
        quality_by_cycle_metrics=... quality_yield_metrics=... picard_{wgs,hs}_metrics=... \
        mosdepth_{global_dist,region_dist,summary}=... mosdepth_{regions,quantized}=... \
        markduplicates_metrics=<dir>/<sg>.markduplicates-metrics \
    --sgs <sg> --meta sequencing_type=<type> ar-guid=<arGuid>
```

where `<dir>` is `[exome/]cram_qc/<sg>`. `cram_qc` is a Metamist analysis type of its own, apart from the
dataset-level `qc` of cpg-flow's MultiQC reports, and has to exist in Metamist before a registering run.

A QC-only SG's `cram_qc` analysis is registered as soon as its QC is done. Its `--output` is the registered
CRAM, from the sample sheet's `reads`. Its duplicate metrics come from the run that aligned it, and are added
only if they're in the SG's dataset at `[exome/]cram_qc/<sg>/` (nf-align) or `[exome/]qc/markduplicates_metrics/`
(cpg-flow). Otherwise the analysis has no `markduplicates_metrics` key, so consumers can't rely on that key;
the other secondary files are always all present.

The sample sheet builder only counts a `cram_qc` analysis whose primary file is the SG's current CRAM, so a
realigned CRAM isn't covered by its predecessor's QC, and a QC-only run never leaves two active `cram_qc`
analyses for one CRAM.

### 8. Failures

- **Per-SG tasks:** an exit code of 137 or 140 (OOM) is retried twice, with memory scaled by the attempt
  number. Any other failure is **ignored**: that SG is dropped and the rest carry on.
- **Registration** is retried once and then ignored.
- **Spot preemption** is retried by `google.batch.maxSpotAttempts = 3`.
- **`BUILD_SAMPLESHEET`**, **`PREFLIGHT`** and the run-records check fail the whole run.
- **No report:** dropped SGs aren't recorded or listed. **The run still exits 0.** Seqera Platform shows the
  failed tasks, and Metamist the analyses the run registered, by `ar-guid`.
- **Recovery is `-resume`:** the sample sheet is cached, so it carries on with the same SGs. It re-runs only
  the failed tasks and those after them, and the pre-flight check is cached, so the files the first attempt
  published don't stop it. A new run (a new session) is different. An SG that failed at fingerprinting or registration has a published CRAM but no analysis, and a new run's
  pre-flight check stops on it until it's registered or deleted by hand. An SG whose QC or QC registration
  failed has a registered CRAM, so a new run takes it as already aligned, and only a `--qc_aligned` run
  retries its QC, once its leftover QC files are deleted.

## Accepted gaps

These follow from the decisions above, and are accepted for now:

- **Concurrent runs:** two runs launched close together for the same SG can both pass the pre-flight
  check, and the second registers the first run's CRAM.
- **One stray file holds up the run:** an unregistered CRAM, or QC files from a failed attempt, stop every
  SG in a new run at the pre-flight check until someone registers or deletes them by hand. Resuming the run
  that left them doesn't stop.
- **No failure report:** a run with dropped SGs exits 0 and records nothing about them. Seqera shows the
  failed tasks.
- **No fingerprint gapfill:** CRAMs registered before nf-align fingerprinted stay without one. Adding one
  needs a process that can update or replace their analyses.
- **QC gapfill is opt-in:** an SG whose QC failed stays without CRAM QC until a `--qc_aligned` run, after
  its leftover QC files are deleted. `--qc_aligned` isn't the default, because every CRAM aligned by
  cpg-flow has QC files but no `cram_qc` analysis, so any cohort run would quietly QC them all.
- **A QC-only analysis can include a file the run didn't write:** the duplicate metrics, found by path.
  This bends the rule that a run only registers files it wrote, for a file that's evidence about the CRAM
  rather than the CRAM itself.
- **QC of a pre-nf-align CRAM** needs a `.dict` beside the DRAGEN-masked FASTA. Without one, Picard fails
  and the SG's QC is dropped.
- **The `cram_qc` analysis type** has to exist in Metamist before a registering run.
- **Alt awareness isn't covered by the stub tests:** their placeholder reference has no `.alt` or `.dict`.
  Only real CRAMs show whether it took effect (RISKS #22).

## Resources (`conf/base.config`)

| process | CPUs | memory | disk (genome / exome) | machine type, by attempt |
|---|---|---|---|---|
| `ALIGN_TO_CRAM`, `ALIGN_SHARD` | 16 | 48 GB × attempt | 400 / 100 GB | `n2-standard-16`, `n2-highmem-16`, `n2-highmem-32` |
| `MERGE_TO_CRAM` | 8 | 24 GB × attempt | 400 / 100 GB | `n2-standard-8`, `n2-highmem-8`, `n2-highmem-16` |
| `SORT_INDEX` | 8 | 32 GB × attempt | 550 / 250 GB | `n2-highmem-8`, `n2-highmem-16`, `n2-highmem-16` |
| `SOMALIER_EXTRACT` | 1 | 4 GB × attempt | 100 / 20 GB | `e2-standard-2`, `e2-standard-4`, `e2-standard-8` |
| `VERIFYBAMID` | 4 | 12 GB × attempt | 100 / 20 GB | `n2-standard-4`, `n2-highmem-4`, `n2-highmem-8` |
| `SAMTOOLS_STATS` | 4 | 4 GB × attempt | 100 / 20 GB | `e2-standard-4`, `e2-standard-4`, `e2-standard-8` |
| `PICARD_METRICS` (two JVMs, 40% heap each) | 4 | 12 GB × attempt | 100 / 20 GB | `n2-standard-4`, `n2-highmem-4`, `n2-highmem-8` |
| `MOSDEPTH` (3 decompression threads) | 4 | 4 GB × attempt | 100 / 20 GB | `e2-standard-4`, `e2-standard-4`, `e2-standard-8` |
| `BUILD_SAMPLESHEET`, `REGISTER_CRAM`, `REGISTER_CRAM_QC` | 1 | 4 / 2 GB | default | `e2-standard-2` |
| `PREFLIGHT` | native `exec:`, in the Nextflow process | | | |

**Every process names one exact `machineType`**, chosen by `task.attempt`, so an OOM retry lands on a
bigger VM. Don't replace these with a list of types or a wildcard (`n2-*`). With
`NXF_CLOUDINFO_ENABLED=false`, Nextflow records a task's machine type in the trace only when it was given
an exact type. Seqera's cost estimates come from that record, so with a list or wildcard Google Batch
still picks a VM, but the type isn't recorded and the task gets no cost estimate ([RISKS.md](RISKS.md)).
Each type leaves memory headroom over the request, because a request equal to the VM's memory may not
fit. The align processes use N2 because bwa-mem3 needs AVX2.

## Files

```
main.nf                       param types and validation, run-records check, sample sheet, per-SG destinations and outputs, pre-flight, ALIGN, somalier, CRAM QC, registration
nextflow.config               param defaults, profiles (standard: Google Batch, docker, test)
nextflow_schema.json          launch-form schema (main.nf does its own validation)
conf/base.config              containers, resources, error strategy, per-SG resource labels, publishDir
subworkflows/local/align/     ALIGN
subworkflows/local/cram_qc/   CRAM_QC
modules/local/<process>/      build_samplesheet, preflight, sort_index, align_to_cram, align_shard, merge_to_cram, somalier_extract, verifybamid, samtools_stats, picard_metrics, mosdepth, register_analysis
bin/build_samplesheet.py      cohorts or SG IDs, from any datasets -> sample sheet (needs only `metamist`)
bin/align_stream.sh           per-kind input -> bwa-mem3 BAM stream
bin/dedup_to_cram.sh          dupblaster -> sort -> CRAM 3.0 + metrics
bin/threads.sh                CPU split between the tools in a pipe
bin/check_alignments.sh       quickcheck + at least one record
test/                         pytest: sample sheet unit tests, stub runs of main.nf
```

## Running

```bash
# production (the launcher adds -work-dir and google.project)
nextflow run main.nf \
    --dataset <analysis ds> --sequencing_type genome --arGuid <guid> --cohorts COH1,COH2 \
    --align_container <image> --cpg_container <image>
  
# or realistically... 
analysis-runner seqera ...

# recovery: the same params, with -resume
nextflow run main.nf -resume ...

# offline wiring check. Stub runs have no Metamist: sheet.config sets
# process.withName:BUILD_SAMPLESHEET.ext.stub_samplesheet to a prepared sheet, or the sheet is empty
nextflow run main.nf -profile test -stub -c sheet.config --dataset ds --sequencing_type genome \
    --cohorts COH1 --ref_fasta refs/ref.fasta

# tests (sample sheet unit tests + stub runs; the stub runs are skipped without nextflow)
uv run pytest
```

`dupblaster_extra_args` is a string that is passed into `ALIGN_TO_CRAM` and `MERGE_TO_CRAM` and appended
verbatim to the dupblaster command, so quote it as you would for bash:

```bash
--dupblaster_extra_args "--read-name-format 'regex:^(?<su>\w+:\d+):(?<tile>\d+):\d+:\d+'"
```

## Status and next steps

**Done**, and stub-tested: everything described above. That is the multi-dataset design (first implemented
from [PLAN.md](PLAN.md)), fingerprinting with one `cram` analysis per CRAM, CRAM QC, `--qc_aligned`, the
cached pre-flight check, dropping the failure report and removing `--samplesheet` and `crams.csv`.
`nextflow lint` is clean, and the 61 tests pass. The stub runs span two datasets and cover every input kind,
the exome layout,
registration paths and projects (including `accessLevel test`), injected failures, the run-records and
pre-flight clashes, `-resume` (after success, after dropped SGs, and after a failed pre-flight check, each
reusing the cached sample sheet),
cohorts and SG ID modes, a sheet row without a dataset, per-SG resource labels,
QC-only SGs and parameter validation.

**Done**, a real test on real data - alignmets are quick and cheap (compared with DRAGMAP). CRAMs generated using 
this new process have been run through the existing variant calling process, and validated with the original version
of the same SG's gVCF. ~98% concordance detected, with the majority of discrepant events resulting from multiple
mapping, or low quality detections. Given the change in tooling, and the change in reference, some discrepancies were
expected, and this is an acceptable margin.
