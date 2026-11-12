# Ambiguous design choices

These are choices where a reasonable person could have gone the other way. They're recorded so that
changing one is a deliberate act rather than a "fix". Part A was agreed during the design interview.
Part B was decided during implementation and **has not been explicitly confirmed**; please review it.

## A. Agreed during design, but could reasonably flip

| choice                                                                                          | the alternative                                                         | why this way                                                                                                             |
|-------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| Failed SGs are dropped and the run exits 0                                                      | `finish` / `terminate`                                                  | a single bad input shouldn't cost a cohort run                                                                           |
| Failed SGs aren't reported                                                                      | `failed.csv` with the step each SG failed at                            | Seqera already lists the failed tasks, the bookkeeping was a large part of `main.nf`, and `-resume` is the recovery path |
| The pre-flight check is a cached task                                                           | a claims file in the work dir                                           | Nextflow's cache already records that the check passed                                                                   |
| Already-aligned SGs stay in the sheet as `aligned`                                              | skip in this current run                                                | QC follows alignment, and can be run even if the alignment wasn't done in this process                                   |
| Single-unit fast path (`ALIGN_TO_CRAM`) kept separate from shard + merge                        | always shard + merge, even for 1 unit                                   | skips a name-sort and a BAM round trip                                                                                   |
| One read group per SG (`ID=SM=<sg>`, `PL:ILLUMINA`)                                             | one read group per lane (`<sg>.<unit>`)                                 | one library per SG for dupblaster and downstream tools; lane information is lost                                         |
| dupblaster sequencing-duplicate detection left on                                               | off by default                                                          | richer metrics; it fails loudly on odd read names, with an override available                                            |
| dupblaster's clamp exit (code 1) treated as failure                                             | parse stderr and accept it                                              | shouldn't happen with short reads; parsing stderr is fragile                                                             |
| QC metrics are registered as one `cram_qc` analysis per CRAM, with the CRAM as the primary file | secondary files on the `cram` analysis; or the existing `qc` type       | QC can't hold up the CRAM, and `qc` already means dataset-level MultiQC reports                                          |
| All of an SG's QC files go in one `[exome/]cram_qc/<sg>/` folder                                | cpg-flow's per-tool `qc/<tool>/` folders                                | one CRAM's QC is in one place; nothing binds nf-align to the old paths yet                                               |
| CRAM QC is all or nothing, and registered only after the CRAM                                   | register whatever succeeded, as soon as it's ready                      | consumers can rely on the set of files, and a `cram_qc` analysis never points at an unregistered CRAM                    |
| Four QC processes, one per image                                                                | one per tool (cpg-flow's four stages, plus mosdepth); one combined task | stages the CRAM four times, not five, without a custom image                                                             |
| Already-aligned SGs without CRAM QC are QC'd only with `--qc_aligned`                           | always gapfill; never                                                   | cpg-flow's CRAMs have no `cram_qc` analysis, so always gapfilling would QC every one of them in any cohort run           |
| `--cram_qc false` still publishes the duplicate metrics to `cram_qc/<sg>/`                      | publish them only with QC                                               | they come from alignment, and have always been published                                                                 |
| No memory watchdog for OOM                                                                      | monitor RSS and kill early                                              | `pipefail` + record checks + OOM retry already cover the header-only failure mode                                        |

## B. Decided during implementation

* **Align tasks request 48 GB, not the 32 GB agreed.** bazam's 16 GB heap, plus bwa-mem3 at ~15 GB,
   dupblaster at up to 5 GB and the sort buffers, comes to about 38 GB. A 16-vCPU standard machine has
   64 GB, so in practice this shouldn't cost more than the 32 GB request would have.
   `MERGE_TO_CRAM` asks for 24 GB (8 CPUs, so it fits an `*-standard-8`).
* **Registration is retried once, then ignored.** `-resume` retries it again.
* **`BUILD_SAMPLESHEET` fails the run if no SG is in scope** (for example, every SG is the other
   sequencing type), rather than doing an empty run. A sheet where every SG is already aligned is still
   a valid run that does nothing.
* **`--bam` (uncompressed BAM), not `--bam=0`.** dupblaster's README uses `--bam=0`, while bwa-mem3's
   own docs give `--bam` for uncompressed output, so the aligner's docs were followed. Check this against
   the pinned binary.
* **`samtools sort -u` (uncompressed) feeds the CRAM writer**, and `samtools merge -u` feeds dupblaster.
   This saves a compress/decompress step on each pipe. It's a change from the port.
* **FASTQ.ora is streamed, not decompressed to disk.** Each `orad` writes to a FIFO in the background,
    and `align_stream.sh` `wait`s on both after bwa-mem3 exits, so a failed decompression fails the task
    (process substitution would hide it). An EXIT trap kills any `orad` still blocked opening its FIFO
    if bwa-mem3 exits first; otherwise it keeps the output pipe open and the task hangs. Decompressing to
    disk, as the port did, needs ~200 GB for a 30x genome and delays alignment until it finishes.
* **Stubs:** `BUILD_SAMPLESHEET`'s stub copies the sheet named by its `ext.stub_samplesheet`, which the stub
    tests set from a config, or else writes a sheet with only a header (nothing to align). And
    `REGISTER_ANALYSIS`'s stub (as `REGISTER_CRAM` and `REGISTER_CRAM_QC`) echoes the `--dry` command it
    would have run.
* **`SORT_INDEX` stages the input reference under `in/`**, next to the reads. This avoids a name
    collision with the alignment reference when the two are the same FASTA.
* **The `test` profile caps resources at 2 CPUs / 4 GB.** The local executor would otherwise refuse
    the 16-CPU / 48 GB requests.
* **The run-records check looks only for `samplesheet.tsv`**, the one run record left.
* **Every process names one exact machine type per attempt**, not a family or wildcard
    (`n2-*,n2d-*,c2d-*,c3d-*`), which was the earlier design. A wildcard gives Google Batch more choice,
    and so possibly cheaper or more available VMs. But Seqera's cost estimates need an exact type
    ([RISKS #32](RISKS.md)), and an exact type per attempt also makes the OOM-retry VM sizes explicit.
* **`PICARD_METRICS` asks for 12 GB, not the 16 GB agreed.** Its first VM, `n2-standard-4`, has
    exactly 16 GB, and a request equal to the VM's memory may not fit (see `conf/base.config`). Each JVM gets
    a 40% heap, about 4.9 GB, which is below the 8 GB cpg-flow gave each Picard job.
    An OOM retries on `n2-highmem-4` with 24 GB. `VERIFYBAMID` is given 12 GB on the same machine types, for
    the same reason.
* **CRAM QC validates its own params.** With `cram_qc` on, the run fails at launch if
    `verifybamid_container`, `picard_container`, `ref_fasta_dict`, or any of the run's sequencing type's
    contamination resources and interval list are unset. Until now `ref_fasta_dict` wasn't checked at all.
* **QC runs even when the SG's fingerprint fails.** QC starts from the aligned CRAM, alongside
    `SOMALIER_EXTRACT`. Such an SG is reported at `somalier`, and its published QC files, like its CRAM, are
    clashes for the next run.
* **`cram_qc` secondary files are ordered as the tools run, then the duplicate metrics.** The order only
    shows on the command line, and in Metamist if it keeps it.
* **An aligned row with an empty `cram_qc` is QC'd** with `--qc_aligned`. An empty column and a missing one
    look the same, but the builder always writes the column, and there's no other way to give a sheet.
* **Param defaults stay in `nextflow.config`; `main.nf` declares only their types.** `conf/base.config`
    reads several params (containers, `arGuid`, `dataset`) while the config is parsed, before `main.nf`'s
    `params {}` block exists, so the defaults can't move into the script. Each param is therefore named in
    both files.
* **`BUILD_SAMPLESHEET` fills in `aligned_ref_fasta`**, so every `aligned` row has a `reads_ref`, and `main.nf`
    never guesses an aligned CRAM's reference.
* **Publishing routes files by name pattern** (`*.cram*` to `cramDir`, `*.markduplicates-metrics` and every QC
    tool's output to `qcDir`), not by a per-file map. A new output file in one of those processes would be
    published without being in `outputs(meta)`, so without a pre-flight check: add it to both.
