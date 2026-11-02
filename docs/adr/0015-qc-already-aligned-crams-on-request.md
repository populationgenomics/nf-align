# CRAM QC of already-aligned SGs, on request

With `--qc_aligned true`, a run also gives **CRAM QC** to each already-aligned SG whose CRAM has no `cram_qc`
analysis. The SG is a **QC-only SG**: it isn't realigned, fingerprinted or re-registered. Only its QC runs, against
the reference its CRAM was written with. A run whose SGs are all aligned does only QC. `--qc_aligned` is off by
default, and needs `--cram_qc` on.

This reverses part of [ADR 0014](0014-cram-qc-one-folder-one-analysis.md), under which an already-aligned SG
never got CRAM QC. That followed the fingerprints in [ADR 0012](0012-one-analysis-per-cram.md), but the two
cases differ. A fingerprint has to go on the existing `cram` analysis, and `cpg_utils.metamist_registration` can
only create analyses. CRAM QC is an analysis of its own, so it can be added to a registered CRAM later.

- **Which SGs**: the sample sheet builder fetches `cram_qc` analyses with the `cram` ones. An aligned row's
  `cram_qc` column is the ID of the latest completed, active `cram_qc` analysis whose primary file
  (`outputs.path`) is that row's CRAM, or empty. Only analyses of the same CRAM count, so a realigned CRAM
  isn't covered by its predecessor's QC. A `--samplesheet` without the column has every aligned SG QC'd.
- **Reference**: the CRAM's own (`reads_ref`, or `--aligned_ref_fasta`), with the `<name>.dict` beside its
  `<name>.fasta`. `ref_fasta`'s is still `--ref_fasta_dict`.
- **Published**: only the QC tools' files, to the same `[exome/]cram_qc/<sg>/`. They go through the pre-flight
  check like every other output, so one left by an earlier, failed attempt stops the run
  ([ADR 0009](0009-never-overwrite-preflight-and-publishdir-no-overwrite.md)).
- **Registered**: one `cram_qc` analysis, with the registered CRAM as its primary file, and no wait for a `cram`
  registration. The duplicate metrics come from alignment, which an earlier run did. If that run left them in
  the SG's dataset, at nf-align's `cram_qc/<sg>/` or at cpg-flow's `qc/markduplicates_metrics/`, they're added
  as the `markduplicates_metrics` secondary file. Otherwise the analysis has no such key.

## Considered Options

- **Always gapfill, with no param**: QC would fill itself in, as alignment does. But CRAMs aligned by cpg-flow
  have QC files and no `cram_qc` analysis, so any cohort run would quietly QC every one of them.
- **QC every aligned SG on request, even those with CRAM QC**: simpler, but would leave two active `cram_qc`
  analyses for one CRAM.
- **Leave the duplicate metrics out of a QC-only analysis**: keeps the run to registering only what it wrote
  ([ADR 0013](0013-only-register-files-this-run-wrote.md)), but leaves out metrics that usually exist.
- **Overwrite, or skip, an SG whose QC files are already there**: either would let a failed QC be retried
  without anyone deleting files by hand. Overwriting goes against ADR 0009, and skipping hides the problem.

## Consequences

- A `cram_qc` analysis may lack `markduplicates_metrics`, so consumers can't rely on that key. Its other
  secondary files are still all or nothing.
- A `cram_qc` analysis can point at a file this run didn't write: the duplicate metrics, found by path. That
  bends ADR 0013's rule, for a file that's evidence about the CRAM, not the CRAM itself.
- An SG that failed at `qc` or `register_qc` can be QC'd again by a `--qc_aligned` run, once its leftover QC
  files are deleted. Before this ADR, nothing could retry it.
- QC of a pre-nf-align CRAM depends on the DRAGEN-masked reference having a `.dict` beside it. If it doesn't,
  Picard fails and the SG is reported at `qc`.
- `failed.csv` reports a QC-only SG at `qc` or `register_qc`. `crams.csv` still lists only the CRAMs the run
  aligned.
