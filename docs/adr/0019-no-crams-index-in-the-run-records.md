# No index of the run's CRAMs: the run records are the sample sheet

`WRITE_RUN_RECORDS` and its `crams.csv` (`sg_id, dataset, cram, crai`, one row per newly aligned SG) are
removed. The run records directory, `alignment_samplesheets/<arGuid>/` in the analysis dataset bucket, now
holds only `samplesheet.tsv`. The run-records check looks for that file alone.

The index repeated what is recorded elsewhere. Each `cram` and `cram_qc` analysis the run registers carries
its `ar-guid` in its meta, so Metamist can list what a run produced. Seqera lists every task, including the
failures. The index was also less exact than it looked. It listed each CRAM's destination as worked out in
`meta`, not a checked publish. It included CRAMs whose fingerprinting or registration then failed, and it left
out QC-only SGs. Keeping it took a run-level process, with its own error strategy and sorting to keep its cache
stable on `-resume`.

This amends the consequence in [ADR 0017](0017-failed-sgs-are-dropped-not-reported.md) that says to compare
`crams.csv` with the sample sheet. Compare the analyses registered under the run's `ar-guid` instead.

## Considered Options

- **Keep `crams.csv`**: one file per run that needs no Metamist query, but it duplicates Metamist and Seqera,
  and isn't exact.
- **Write it only with `--register false`**, where Metamist has nothing: a second code path for test runs,
  which can list the bucket instead.

## Consequences

- To see what a run produced, query Metamist for analyses with its `ar-guid`, or look at the run in Seqera.
- A CRAM that was published but never registered (its fingerprinting or registration failed) is in no
  record. Seqera shows the failed task, and a new run's pre-flight check finds the CRAM
  ([ADR 0013](0013-only-register-files-this-run-wrote.md)).
- A `--register false` run leaves no record of its CRAMs, only the files in the buckets.
