# Risks: things that may come back to bite

Each entry gives the risk, why it matters, and what currently guards against it. A ✅ marks a risk that
was verified by experiment or by reading source during design; anything without a ✅ is inferred from
documentation.

## Nextflow behaviour

* **`publishDir overwrite: false` skips silently.** ✅ (Nextflow 26.10.4, local) An existing target is
   left in place with no warning, and the exit code is 0 (`workflow.output.overwrite = false` behaved the
   same on 26.04). Overwriting isn't allowed in our buckets, so `false` is required, and the pre-flight
   check is the only thing that catches a clash. Don't remove it.
* **Publishing is asynchronous.** ✅ (with workflow outputs; assumed to be the same for `publishDir`)
   A downstream task can start before its input has been published, and locally a partly copied file
   was observed. Registration polls because of this. A cross-bucket copy of a large WGS CRAM
   could take longer than `register_wait_minutes`. If that happens, the SG's registration fails while its
   CRAM is in place. `-resume` registers it; a new run trips the pre-flight check on it.
* **`Path.toString()` drops the bucket on GCS.** ✅ `file('gs://b/x/y').toString()` gives `/x/y`.
   `uriOf()` in `main.nf` uses `toUriString()` for cloud paths. Any new code that builds a `gs://`
   string from a Path must go through `uriOf()`.
* **`error()` inside a channel operator is hidden.** ✅ The console shows only "Unexpected error
   [InvocationTargetException]". Use `fail()` in `main.nf` (log, then error) for any check that runs
   inside an operator.
* **`--param ''` on the command line arrives as boolean `true`, not empty.** ✅ A param typed `String`
   in `main.nf` gets the string `'true'`, which the required-param check accepts. `validateParams()` catches
   it only for `dupblaster_extra_args`.
* **A failed run exits 0, and says nothing about it.** ✅ Under `ignore`, see
   [DESIGN.md](DESIGN.md#8-failures). Launchers, Seqera Platform and any scripts
   watching the exit code will report success, and there's no `failed.csv`. Look at the failed tasks in
   Seqera, or compare the analyses registered with this `ar-guid` in Metamist with the sample sheet.
* **Typed params need Nextflow 26.04.** ✅ (25.10.0, 26.04.4, 26.10.4, local) The `params {}` block in
   `main.nf` uses `String?` for optional params, which 25.10 can't parse, so `nextflowVersion` is `>=26.04`.
   Config values aren't type-checked, only command-line ones: `cohorts` declared `String?` still takes a list
   from a config file.
* **Shard count is part of the task hash.** Changing `realignment_shards` between `-resume` runs
    re-aligns every CRAM/BAM input from scratch.
    - **So are each SG's destinations.** `meta` carries the SG's Metamist project and output folders, and
      `meta` is part of every per-SG task hash, and of `PREFLIGHT`'s. Changing `output_root`, `accessLevel`
      or `sequencing_type` between `-resume` attempts re-runs alignment for every SG.

## Re-running after a partial failure

* **A new run rebuilds the sheet; a resumed one doesn't.** After a run where some SGs failed:
    - **A new run** (new `arGuid`): the successful SGs now have CRAM analyses, so they come back as `aligned`,
      and only the failures are aligned again.
    - **Registration failures:** if registration failed, the CRAM and its fingerprint exist but aren't
      registered, so a new run trips the pre-flight check on them until they're registered or deleted by hand.
    - **Fingerprint failures:** likewise, the CRAM is published but not registered, so the next run trips the
      pre-flight check on it.
    - **`-resume`** is the way to recover: a passed pre-flight check is cached, so the earlier attempt's
      outputs aren't clashes, and `publishDir` skips them. Failed tasks, including registrations,
      run again. The cached sample sheet is reused, so SGs registered by the earlier attempt are republished
      from cache and not registered again. If the sheet's inputs change between attempts (cohorts, SG IDs,
      `aligned_ref_fasta`), it's rebuilt from Metamist, and the pre-flight check runs again and fails. The cache lives
      in the work dir, so cleaning the work dir (or resuming from a different one) brings the check back, and
      it then fails on the earlier attempt's files.

## Tools

* **bwa-mem3 is young and moving fast.** v0.15.0 was released on 2026-10-14, with four releases in a
    month, and it was forked from bwa-mem2 and renamed in April 2026. Pin the version in the image, and
    recheck flags (`--bam`, `-K`, `-p`) whenever the pin moves. The design relies on documented
    behaviour that hasn't been run here:
    - interleaved FASTQ on stdin (`-p -`)
    - output that is query-grouped (`GO:query`)
    - `-K` making output identical regardless of thread count
* **bwa-mem3 needs AVX2.** On a CPU without it, it exits with code 2. N1 machines may not have AVX2,
    which is why the align processes use N2 machine types. If you change their machine family, check
    that it has AVX2.
* **An OOM-killed aligner can leave a valid but empty BAM** that exits 0 when its exit status is lost
    (no `pipefail`, or a background producer whose status nobody collects). This design uses `pipefail`,
    pipes bazam into bwa-mem3, and runs `check_alignments.sh`. FASTQ.ora is the one FIFO path: each
    background `orad` is `wait`ed on explicitly, and an EXIT trap kills them if bwa-mem3 exits first.
    Any new FIFO or background process needs both.
* **dupblaster exits 1 for a "coordinate clamped" warning**, with complete output. ✅ (v0.3.0 source,
    `src/main.rs:804`) The exit code is the same as for real errors, and only stderr tells them apart
    (`(Exiting non-zero.)` vs `Premature exit`). It shouldn't happen with ~150 bp reads under the default
    `--max-read-length 1000`. If it does, the SG fails; rerun it with a larger `--max-read-length` through
    `dupblaster_extra_args`.
* **dupblaster fails on read names that aren't Illumina-style.** ✅ This is because sequencing-duplicate
    detection is on by default. SRA-renamed, MGI or simulated reads, and some re-headed CRAMs, will fail.
    Override with `--sequencing-duplicate-detection off` or `--read-name-format regex:...`.
* **dupblaster only checks ordering heuristically.** ✅ Input that isn't query-grouped usually errors,
    but single-end data with no secondary or supplementary records gives silently wrong calls. Anything
    new feeding dupblaster must stay query-grouped: aligner output, `samtools sort -n`, or `merge -n`.
* **dupblaster's temp disk grows with depth.** ✅ About 5 GB at 30x and about 50 GB at 300x, under
    `--tmp-dir tmp` in the task directory. High-depth exomes or genomes may need more disk.
* **bazam's heap is capped at 16 GB** (`-Xmx16g`). Together with bwa-mem3 (~15 GB at 16 threads),
    dupblaster (1–5 GB) and the sort buffers, that's why align tasks request 48 GB
    ([AMBIGUOUS-CHOICES.md](AMBIGUOUS-CHOICES.md)).

## Reference and index

* **bwa-mem3 index sidecars.** If `<prefix>.hdr` or `<baseprefix>.dict` sits next to the index, bwa-mem3
    copies its `@SQ` lines verbatim. `ref_fasta_dict` is staged deliberately, and is the alt-annotated
    `.dict` (`samtools dict --alt`), so the alt contigs keep `AH:*`. Replacing it with a plain
    GATK `.dict` would drop those tags from every CRAM header.
* **How alt awareness affects alignment wasn't checked.** We align alt-aware (`<ref_fasta>.alt` is staged),
    so reads in alt-haplotype regions should stay primary on the primary assembly instead of
    getting MAPQ 0. The stub tests can't show this. Check the first real CRAMs: their `@SQ` alt contigs
    carry `AH`, and MAPQ in alt-haplotype regions isn't pulled down to 0.
* **The input reference falls back to the alignment reference.** ✅ A CRAM whose assay has no
    `reference_assembly.location` used to be decoded with `ref_fasta`, and failed in bazam with a reference
    MD5 mismatch when it was written against another reference. Such a CRAM now fails the sample sheet build,
    and a `cram` row with no `reads_ref` fails the run before any work starts.

## Metamist and GCS

* **`cpg_utils.metamist_registration` only accepts `gs://` paths that exist.** ✅ Registration can't
    run locally, which is why the `test` profile sets `register = false`.
* **Two runs for the same SG.** Two runs launched close together for the same SG can both pass the
    pre-flight check. The second run's publish is skipped silently, and it registers the first run's CRAM
    under its own `ar-guid`. That leaves two active `cram` analyses, and the sample sheet builder takes
    the latest. This is accepted.
* **Per-SG resource labels only move cost under `seqr`.** ✅ A `resourceLabels` closure is evaluated
    per task against `meta` (Nextflow 26.10.4, local executor), and
    `test_per_sg_tasks_are_labelled_with_the_sgs_dataset` checks every per-SG task's label in a stub run.
    Production billing for a run under `seqr` shows Google Batch applies the labels, and each SG's
    compute appears under its own dataset. Only the seqr GCP project's costs are redistributed by the
    `dataset` label, though. A run under any other analysis dataset is billed wholly to that dataset,
    and the label only tags the rows (DESIGN.md §4).

## Cost

* **Seqera's cost estimates need exact machine types.** With `NXF_CLOUDINFO_ENABLED=false`, Nextflow
    only records a task's machine type in the trace when the process names one exact `machineType`.
    For a list of types or a wildcard, Google Batch picks the VM, but the type isn't recorded, and Seqera
    has nothing to estimate the task's cost from. Nothing fails, so the estimates just go missing. That's
    why every process in `conf/base.config` names one exact type per attempt. Keep it that way when
    adding a process or changing resources. This hasn't been checked against a Seqera run yet.

## CRAM QC

* **Failed CRAM QC is retried only by hand.** An SG that fails at `qc` or `register_qc` has a registered
    CRAM, so no later run aligns it again, and its published QC files are left unregistered. A `--qc_aligned`
    run QCs it again, but its pre-flight check stops on those files until someone deletes them.
* **Path-based consumers won't find the QC files.** cpg-flow's `CramMultiQC` and anything else that
    reads `qc/<tool>/` paths need to read the `cram_qc` analysis or `cram_qc/<sg>/` instead.
* **A QC-only SG's CRAM is QC'd against the DRAGEN-masked reference's `.dict`.** ✅ `Homo_sapiens_assembly38_masked.dict`
    sits beside `aligned_ref_fasta`. It has the same 3,366 contigs and lengths as `ref_fasta`'s `.dict`, but
    different md5s on 786 of them, all alt contigs. CollectHsMetrics compares the CRAM's header with the
    target and bait lists' headers, md5s included (Picard 3.5.0, `CollectTargetedMetrics`). The regions are
    stored as BED, with no header, and `PICARD_METRICS` makes the interval list from the CRAM's own `.dict`
    for both sequencing types, so neither metric can see a mismatched header. Every reference a CRAM is QC'd
    against needs a `.dict` beside it.
