# CRAM QC: one folder and one `cram_qc` analysis per CRAM

Each newly aligned CRAM gets **CRAM QC**: VerifyBamID, `samtools stats`, Picard CollectMultipleMetrics,
and Picard CollectWgsMetrics (genome) or CollectHsMetrics (exome). These are cpg-flow's `CramQc*` stages,
run in the `CRAM_QC` subworkflow. All of an SG's metrics files, including the duplicate metrics from
alignment, are published to one folder, `[exome/]cram_qc/<sg>/`. They're registered as one Metamist
`cram_qc` analysis, with the CRAM as the primary file and every metrics file as a secondary file.

cpg-flow wrote each tool's output to its own folder (`qc/verify_bamid/`, `qc/samtools_stats/`, five
`qc/<metric>/` folders for CollectMultipleMetrics, `qc/picard_{wgs,hs}_metrics/`). That layout was a legacy
of production-pipelines and isn't followed here. nf-align hasn't published at scale, so no existing outputs
tie it to the old paths. The file names inside the folder are still cpg-flow's (for example
`<sg>.verify-bamid.selfSM`), so MultiQC's filename patterns keep matching.

CRAM QC is kept apart from the CRAM's own registration. A QC failure never stops the `cram` analysis from
being registered. The `cram_qc` analysis is registered only after the `cram` analysis has been, so it can't
point at an unregistered CRAM. It's all or nothing: if any QC task fails, no `cram_qc` analysis is
registered.

`cram_qc` is a new Metamist analysis type. It isn't cpg-flow's `qc`, which records dataset-level MultiQC
reports.

## Considered Options

- **cpg-flow's per-tool folders**: existing path-based consumers would keep working, but one CRAM's QC
  would be spread over nine folders.
- **Analysis type `qc`, with a meta field to tell it apart**: no new type, but every `type=qc` query would
  have to filter out the dataset-level MultiQC records.
- **Register partial QC**: some metrics would survive one tool's failure, but consumers couldn't rely on
  the set of secondary files.
- **QC files as secondary files on the `cram` analysis**: one record, but QC would gate the CRAM's
  registration, and a QC failure would leave the CRAM unregistered (ADR 0013).

## Consequences

- Path-based consumers written for cpg-flow, such as its `CramMultiQC` stage, can't find these outputs
  where they expect them. They need to read the `cram_qc` analysis, or the new folder.
- The duplicate metrics move from `qc/markduplicates_metrics/` to `cram_qc/<sg>/`, and are published by
  alignment before QC runs.
- An SG that fails at `qc` or `register_qc` has a registered CRAM but leaves published, unregistered
  QC files. The CRAM is registered, so later runs treat the SG as already aligned. They don't realign it,
  don't trip on its QC files, and don't retry its QC, so it stays without CRAM QC until something
  gapfills it.
- Already-aligned SGs don't get CRAM QC, as they don't get fingerprints (ADR 0012).
- The `cram_qc` analysis type has to exist in Metamist before a registering run.
