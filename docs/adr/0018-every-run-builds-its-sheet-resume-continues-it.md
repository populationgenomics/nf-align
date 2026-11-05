# Every new run builds its sample sheet; a resumed run carries on with the one it built

`--samplesheet`, which took a prepared TSV and skipped the Metamist query, is removed. A new run takes
`--cohorts` or `--sgids`, and `BUILD_SAMPLESHEET` builds the sheet and publishes it to the run records. A run
that didn't finish is recovered with `-resume`. `BUILD_SAMPLESHEET` is cached like any other task, so the
resumed run reuses the sheet its first attempt built. It doesn't query Metamist again, where the SGs the first
attempt registered would now show up as aligned. With the same rows, `PREFLIGHT` is cached too
([ADR 0016](0016-preflight-check-is-a-cached-task.md)), and only the tasks that failed or never ran, and those
after them, run.

A prepared sheet was a second way in that the builder's rules didn't cover. It could leave out `dataset`,
`reads_ref` or `cram_qc`, and `main.nf` had to check for that or assume a default. Its run records had no
`samplesheet.tsv`. And it didn't recover by itself: rerunning the same sheet after a partial failure clashed
in the pre-flight check on every SG that had succeeded. Its uses were stub runs, tiny test runs and re-running
a known sheet. `-resume` covers the last of these. For stub runs, a test config can give `BUILD_SAMPLESHEET`'s
stub a prepared sheet (`ext.stub_samplesheet`).

## Considered Options

- **Keep `--samplesheet`**: a way around Metamist, but a second set of input rules to keep in step with the
  builder, and no recovery story of its own.
- **On `-resume`, read the published `samplesheet.tsv` from the run records instead of the cached task**:
  this would survive a lost work dir. But the cached pre-flight check is in the work dir too, so the run would
  fail anyway. And if the cohorts changed between attempts, it would ignore the change silently.

## Consequences

- Every run queries Metamist, so a run needs Metamist access and in-scope SGs, even for a small test.
- Every run publishes `samplesheet.tsv`, so the run-records check always has it to find.
- A resumed run sees the SGs as they were when its first attempt started. An SG aligned somewhere else in
  the meantime is aligned again. Launch a new run (a new `arGuid`) to pick up Metamist's current state.
- `BUILD_SAMPLESHEET`'s cache entry depends on the access level, sequencing type, cohorts or SG IDs,
  `aligned_ref_fasta` and the run-records directory. If one of these changes between attempts, the sheet is
  rebuilt from Metamist. The new rows then miss `PREFLIGHT`'s cache, and the check fails on the earlier
  attempt's files ([RISKS #10](../RISKS.md)).
