# Risks: things that may come back to bite

Each entry gives the risk, why it matters, and what currently guards against it. A ✅ marks a risk that
was verified by experiment or by reading source during design; anything without a ✅ is inferred from
documentation.

## Nextflow behaviour

1. **`workflow.output.overwrite = false` skips silently.** ✅ (Nextflow 26.04) There's no warning and
   the exit code is 0, while the "Outputs" summary still lists the file. That's why the pre-flight
   check exists ([ADR 0002](adr/0002-preflight-check-instead-of-overwrite-false.md)). Don't "simplify"
   this back to `overwrite = false`.
2. **Publishing is asynchronous.** ✅ A downstream task can start before its input has been published,
   and locally a partly copied file was observed. Registration polls because of this
   ([ADR 0005](adr/0005-registration-as-a-polling-process.md)). A cross-bucket copy of a large WGS CRAM
   could take longer than `register_wait_minutes`. If that happens, the SG is reported as failed at
   `register` while its CRAM is in place, and the next run then trips the pre-flight check on it.
3. **A failed index-file write does not fail the run.** ✅ When `crams.tsv` was used (index files must
   be `.csv`, `.json` or `.yaml`), Nextflow logged `Failed to invoke workflow.onComplete event handler`
   and still reported SUCCESS. The stub tests assert that the index files exist.
4. **`Path.toString()` drops the bucket on GCS.** ✅ `file('gs://b/x/y').toString()` gives `/x/y`.
   `uriOf()` in `main.nf` uses `toUriString()` for cloud paths. Any new code that builds a `gs://`
   string from a Path must go through `uriOf()`.
5. **`error()` inside a channel operator is hidden.** ✅ The console shows only "Unexpected error
   [InvocationTargetException]". Use `fail()` in `main.nf` (log, then error) for any check that runs
   inside an operator.
6. **`--param ''` on the command line arrives as boolean `true`, not empty.** ✅ The required-param
   check would accept it.
7. **A failed run exits 0.** ✅ Under `ignore`, see [ADR 0007](adr/0007-ignore-failed-sgs-and-report.md).
   Launchers, Seqera Platform and any scripts watching the exit code will report success.
8. **`failed.csv` is only written if the run reaches the end.** If the sample sheet step or the
   pre-flight check fails, or the run is cancelled, there's no failed-SG report, and the run itself
   fails anyway.
9. **The run-records path is keyed on `arGuid`, and `overwrite = true`.** Two runs with the same
   `arGuid` overwrite each other's `samplesheet.tsv`, `crams.csv` and `failed.csv`. The pre-flight check
   doesn't cover run records.
10. **Shard count is part of the task hash.** Changing `realignment_shards` between `-resume` runs
    re-aligns every CRAM/BAM input from scratch.

## Re-running after a partial failure

11. **Cohorts and SG ID modes recover by themselves; `--samplesheet` mode doesn't.** After a run where some SGs failed:
    - **Cohorts / SG ID mode:** the successful SGs now have CRAM analyses, so they come back as `aligned`, and
      only the failures are aligned again.
    - **`--samplesheet` mode:** the successful SGs' CRAMs now exist, so the pre-flight check fails the
      rerun. Remove those rows first, or rebuild the sheet from cohorts or SG IDs.
    - **Registration failures:** if registration failed, the CRAM exists but isn't registered, so either
      mode trips the pre-flight check until it's registered by hand.
    - **`-resume`** works in both modes: outputs claimed by the session's earlier attempt are allowed to
      exist (ADR 0002). In cohorts and SG ID modes the cached sample sheet is reused, so SGs registered by the
      earlier attempt are republished from cache and not registered again. The claims file lives in the
      work dir, so cleaning the work dir (or resuming from a different one) brings the check back.

## Tools

12. **bwa-mem3 is young and moving fast.** v0.15.0 was released on 2026-10-14, with four releases in a
    month, and it was forked from bwa-mem2 and renamed in April 2026. Pin the version in the image, and
    recheck flags (`--bam`, `-K`, `-p`) whenever the pin moves. The design relies on documented
    behaviour that hasn't been run here:
    - interleaved FASTQ on stdin (`-p -`)
    - output that is query-grouped (`GO:query`)
    - `-K` making output identical regardless of thread count
13. **bwa-mem3 needs AVX2.** On a CPU without it, it exits with code 2. N1 machines may not have AVX2,
    which is why the align processes are limited to `n2-*,n2d-*,c2d-*,c3d-*`. If you add a machine
    family, check it.
14. **An OOM-killed aligner can leave a valid but empty BAM** that exits 0 when its exit status is lost
    (no `pipefail`, or a background producer whose status nobody collects). This design uses `pipefail`,
    pipes bazam into bwa-mem3, and runs `check_alignments.sh`. FASTQ.ora is the one FIFO path: each
    background `orad` is `wait`ed on explicitly, and an EXIT trap kills them if bwa-mem3 exits first.
    Any new FIFO or background process needs both.
15. **dupblaster exits 1 for a "coordinate clamped" warning**, with complete output. ✅ (v0.3.0 source,
    `src/main.rs:804`) The exit code is the same as for real errors, and only stderr tells them apart
    (`(Exiting non-zero.)` vs `Premature exit`). It shouldn't happen with ~150 bp reads under the default
    `--max-read-length 1000`. If it does, the SG fails; rerun it with a larger `--max-read-length` through
    `dupblaster_extra_args`.
16. **dupblaster fails on read names that aren't Illumina-style.** ✅ This is because sequencing-duplicate
    detection is on by default. SRA-renamed, MGI or simulated reads, and some re-headed CRAMs, will fail.
    Override with `--sequencing-duplicate-detection off` or `--read-name-format regex:...`.
17. **dupblaster only checks ordering heuristically.** ✅ Input that isn't query-grouped usually errors,
    but single-end data with no secondary or supplementary records gives silently wrong calls. Anything
    new feeding dupblaster must stay query-grouped: aligner output, `samtools sort -n`, or `merge -n`.
18. **dupblaster's temp disk grows with depth.** ✅ About 5 GB at 30x and about 50 GB at 300x, under
    `--tmp-dir tmp` in the task directory. High-depth exomes or genomes may need more disk.
19. **bazam's heap is capped at 16 GB** (`-Xmx16g`). Together with bwa-mem3 (~15 GB at 16 threads),
    dupblaster (1–5 GB) and the sort buffers, that's why align tasks request 48 GB
    ([AMBIGUOUS-CHOICES.md](AMBIGUOUS-CHOICES.md)).

## Reference and index

20. **No public bwa-mem3 index exists for this reference.** ✅ The Broad bucket only has bwa 0.7 indexes
    and GENCODE bwa-mem2 builds. Building one needs about 76 GB of RAM: run `bwa-mem3 index` on a
    highmem machine once, and store the result next to the FASTA.
21. **bwa-mem3 index sidecars.** If `<prefix>.dict` or `<prefix>.hdr` sits next to the index, bwa-mem3
    copies its `@SQ` lines verbatim. Only the suffixes in `ref_suffixes` (main.nf) are staged
    (`fai, bwt.2bit.64, amb, ann, pac`), so no sidecar is seen. Adding `dict` there would change the CRAM headers.
22. **How the masked FASTA treats alt contigs wasn't checked.** We align without alt awareness (no
    `.alt`). Check how alt contigs show up in the first real CRAMs.
23. **The input reference falls back to the alignment reference.** A CRAM whose assay has no
    `reference_assembly.location` is decoded with `ref_fasta`. If it was actually encoded against
    another reference, decoding fails in `SORT_INDEX` or in bazam.

## Metamist and GCS

24. **Metamist access from Google Batch is untested.** Both `BUILD_SAMPLESHEET` and `REGISTER_ANALYSIS`
    need the task service account to reach Metamist.
25. **`cpg_utils.metamist_registration` only accepts `gs://` paths that exist.** ✅ Registration can't
    run locally, which is why the `test` profile sets `register = false`.
26. **Already-aligned rows are never checked.** If a registered CRAM has been deleted, nothing notices
    until something (later, genotyping) tries to read it.
27. **An analysis with `output` set but its file missing** is treated as already aligned, as in (26).
