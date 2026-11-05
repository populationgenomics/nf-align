# One analysis per CRAM, with the fingerprint as a secondary file

Each newly aligned CRAM is registered as a single Metamist `cram` analysis. Its `.crai` (`index`) and
`.somalier` (`somalier`) are secondary files, and the reference it was written against is `ref_fasta` in its
meta. There is no `somalier` analysis type. `REGISTER_ANALYSIS` runs once per SG, after `SOMALIER_EXTRACT`,
and waits for all three files to be published.

So a CRAM whose fingerprint failed isn't registered. It's reported at `somalier`, and its published CRAM
stops the next run at the pre-flight check ([ADR 0013](0013-only-register-files-this-run-wrote.md)).

Already-aligned SGs are no longer fingerprinted (**gapfilled**). Their CRAM analysis already exists, and
`cpg_utils.metamist_registration` can only create analyses, not add a secondary file to one. The sample sheet
still records an aligned SG's fingerprint, read from its CRAM analysis's `somalier` secondary file, but
`main.nf` doesn't use it.

This replaces the earlier design on this branch, where a CRAM and its fingerprint were registered as separate
`cram` and `somalier` analyses, in their own tasks.

## Considered Options

- **Separate `cram` and `somalier` analyses**: a CRAM was registered without waiting for its fingerprint, but
  a consumer had to join the two records, and `failed.csv` reported up to two rows per SG.
- **Gapfill with a fresh `cram` analysis**: register the existing CRAM again, with the new fingerprint as a
  secondary file. That leaves two active analyses for one CRAM.
- **Gapfill by publishing only**: the next run couldn't see the fingerprint, would try again, and would trip
  the pre-flight check on the `.somalier` it left.

## Consequences

- CRAMs registered before nf-align fingerprinted have no fingerprint, and this pipeline won't make one.
  Fingerprinting them needs a separate process that can update or replace their analyses.
- A CRAM is registered later than before, after its fingerprint, and a fingerprint failure leaves it
  unregistered for someone to sort out by hand.
- `failed.csv` has one row per SG: `align`, `somalier` or `register`.
- The sample sheet builder only queries `cram` analyses, and expects Metamist to return a secondary file at
  `outputs.secondary_files.<type>.path`.
