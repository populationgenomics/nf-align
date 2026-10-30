# Ambiguous design choices

These are choices where a reasonable person could have gone the other way. They're recorded so that
changing one is a deliberate act rather than a "fix". Part A was agreed during the design interview.
Part B was decided during implementation and **has not been explicitly confirmed**; please review it.

## A. Agreed during design, but could reasonably flip

| choice | the alternative | why this way |
|---|---|---|
| Hash shards, not region shards ([ADR 0003](adr/0003-read-name-hash-shards-not-regions.md)) | region shards, as the original brief said ("roughly by regions") | keeps mates together and covers unmapped reads |
| Failed SGs are dropped and the run exits 0 ([ADR 0007](adr/0007-ignore-failed-sgs-and-report.md)) | `finish` / `terminate` | a single bad input shouldn't cost a cohort run |
| Already-aligned SGs stay in the sheet as `aligned` | leave them out until genotyping exists | genotyping is next, and they add no tasks |
| Single-unit fast path (`ALIGN_TO_CRAM`) kept separate from shard + merge | always shard + merge, even for 1 unit | skips a name-sort and a BAM round trip |
| One read group per SG (`ID=SM=<sg>`, `PL:ILLUMINA`) | one read group per lane (`<sg>.<unit>`) | one library per SG for dupblaster and downstream tools; lane information is lost |
| dupblaster sequencing-duplicate detection left on | off by default | richer metrics; it fails loudly on odd read names, with an override available |
| dupblaster's clamp exit (code 1) treated as failure | parse stderr and accept it | shouldn't happen with short reads; parsing stderr is fragile |
| Metrics are published but not registered | register them as their own analysis type | nothing reads them from Metamist yet |
| Each SG's bucket derived from its dataset + `accessLevel`; `output_root` only for local/test runs | the launcher passes an output location | the bucket can't disagree with the SG or the access level |
| An SG's dataset is its sample's project ([ADR 0008](adr/0008-sgs-from-several-datasets-published-with-publishdir.md)) | the project of the cohort that brought it in | the data belongs with its sample, whoever curated the cohort |
| No allow-list of datasets a run may write to | a `--datasets` allow-list | the SG's project picks the bucket, so a wrong cohort can't send data to the wrong place |
| An SG at the wrong access level fails the sample sheet step | skip it with a warning | a test run must never write to a `-main` bucket, and a full run never to a `-test` one |
| Still one sequencing type per run | per-SG `exome/` prefix, mixed runs | resources, shard defaults and genotyping depend on it; a separate change if ever wanted |
| `publishDir` for everything, including run records | keep `output {}` for run records | one publishing mechanism, not two |
| `-resume` skips the run-records check | record run records in the claims file | simpler, and a resume keeps its `arGuid` by definition |
| The concurrent-run gap is accepted ([ADR 0009](adr/0009-never-overwrite-preflight-and-publishdir-no-overwrite.md)) | re-check existence just before each publish | a small race would remain anyway, and it adds a moving part |
| No memory watchdog for OOM | monitor RSS and kill early | `pipefail` + record checks + OOM retry already cover the header-only failure mode |
| bwa-mem3 index built outside the pipeline | an optional `BWA_MEM3_INDEX` process | needs ~76 GB once; not per-run work |

## B. Decided during implementation: please review

1. **Align tasks request 48 GB, not the 32 GB agreed.** bazam's 16 GB heap, plus bwa-mem3 at ~15 GB,
   dupblaster at up to 5 GB and the sort buffers, comes to about 38 GB. A 16-vCPU standard machine has
   64 GB, so in practice this shouldn't cost more than the 32 GB request would have.
   `MERGE_TO_CRAM` asks for 24 GB (8 CPUs, so it fits an `*-standard-8`).
2. **The index files are CSV, not TSV.** Originally this was because Nextflow only accepts `.csv`,
   `.json` or `.yaml` for workflow output index files. `WRITE_RUN_RECORDS` could now write TSV, but the
   names `crams.csv` and `failed.csv` are kept so anything that already reads them keeps working. The
   sample sheet itself is still `samplesheet.tsv`.
3. **`failed.csv` has a `failed_at` column** (`align` or `register`), so a failed registration isn't
   confused with a failed alignment.
4. **Registration is retried once, then ignored** (and reported in `failed.csv`).
5. **`BUILD_SAMPLESHEET` fails the run if no SG is in scope** (for example, every SG is the other
   sequencing type), rather than doing an empty run. A sheet where every SG is already aligned is still
   a valid run that does nothing.
6. **The sheet builder takes `--access-level`** (agreed in the multi-dataset design). It rejects SGs
   whose project doesn't match it, and writes each SG's `dataset` without `-test`. `main.nf` works out
   each SG's Metamist project for registration from `meta.dataset` and `accessLevel`. This replaces the
   earlier design, where `main.nf` passed a single project name and the script had no `accessLevel`
   logic.
7. **Registration `--meta` has `sequencing_type` and `ar-guid` only.** A `source=nf-align` key was
   suggested early in the interview but wasn't in the final summary, so it's left out. It's one line to
   add in `modules/local/register_analysis/main.nf`.
8. **`--bam` (uncompressed BAM), not `--bam=0`.** dupblaster's README uses `--bam=0`, while bwa-mem3's
   own docs give `--bam` for uncompressed output, so the aligner's docs were followed. Check this against
   the pinned binary.
9. **`samtools sort -u` (uncompressed) feeds the CRAM writer**, and `samtools merge -u` feeds dupblaster.
   This saves a compress/decompress step on each pipe. It's a change from the port.
10. **FASTQ.ora is streamed, not decompressed to disk.** Each `orad` writes to a FIFO in the background,
    and `align_stream.sh` `wait`s on both after bwa-mem3 exits, so a failed decompression fails the task
    (process substitution would hide it). An EXIT trap kills any `orad` still blocked opening its FIFO
    if bwa-mem3 exits first; otherwise it keeps the output pipe open and the task hangs. Decompressing to
    disk, as the port did, needs ~200 GB for a 30x genome and delays alignment until it finishes.
11. **The pre-flight check also covers `.crai` and metrics**, not just the CRAM. A half-published SG
    (index present, CRAM missing) still counts as a clash.
12. **Run records live under `alignment_samplesheets/<arGuid>/` in the analysis dataset's bucket**, as
    you specified for the sample sheet, so `crams.csv` and `failed.csv` sit next to it.
13. **The `gcb` container check only applies when the profile list includes `gcb`.** The `docker`
    profile without images would fail at task launch instead.
14. **`dupblaster_extra_args` is one string, inserted verbatim into the dupblaster command line.** Bash
    splits and unquotes it, so it can hold several arguments, but quoting is the caller's job. A bare
    `--dupblaster_extra_args` (which arrives as `'true'`) fails parameter validation.
15. **`cohorts` and `sgids` take commas or spaces** (`COH1,COH2`), or a list from config.
    - **`sgids` fails on an SG that isn't short-read `sequencing_type`**, where cohorts mode skips it with a
      warning. A named SG is an explicit request, so a silent skip would leave it unaligned unnoticed.
16. **Stubs:** `BUILD_SAMPLESHEET`'s stub writes a sheet with only a header (nothing to align), and
    `REGISTER_ANALYSIS`'s stub echoes the `--dry` command it would have run.
17. **`SORT_INDEX` stages the input reference under `in/`**, next to the reads. This avoids a name
    collision with the alignment reference when the two are the same FASTA.
18. **The `test` profile caps resources at 2 CPUs / 4 GB.** The local executor would otherwise refuse
    the 16-CPU / 48 GB requests.
19. **`crams.csv` lists the destination URIs of newly aligned CRAMs**, as worked out in `meta`. It doesn't
    check that the publish succeeded, or that it wasn't skipped. This matches what the workflow output
    index did before.
20. **The run-records check fails if any one of the three files exists**, even when the other two are
    missing. In `--samplesheet` mode `samplesheet.tsv` is never written, but the check still covers the
    other two.
21. **`WRITE_RUN_RECORDS` is a native `exec:` process**, not a `script:` process on the local executor.
    A spike on 26.10.4 showed that a file written to `task.workDir` from Groovy, and declared as a `path`
    output, is published by `publishDir`. It needs no container and no VM.
22. **`WRITE_RUN_RECORDS` fails the run** (`errorStrategy 'terminate'`), as `BUILD_SAMPLESHEET` does.
    Under the default `ignore`, a failure would have lost the run records without notice.
23. **`crams.csv` and `failed.csv` rows are sorted by SG ID.** The records reach the process in task
    completion order, so without sorting the input would change between attempts, and `-resume` would
    rerun the task.
24. **The access-level check comes before the sequencing-type filter.** An SG at the wrong access level
    fails the sample sheet step, even if it would otherwise have been skipped as the other sequencing
    type or as long-read.
25. **`BUILD_SAMPLESHEET` is tagged with the access level**, now that it isn't given a project.
