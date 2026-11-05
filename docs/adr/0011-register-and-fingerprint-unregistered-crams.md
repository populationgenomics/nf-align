# An unregistered CRAM is fingerprinted and registered, not a clash

**Status: superseded by [ADR 0013](0013-only-register-files-this-run-wrote.md).** An unregistered CRAM is a
clash again, and nf-align only registers files it wrote.

An earlier run can publish an SG's CRAM and `.crai` but never register them, for example because the run
failed at registration, or predates fingerprinting. Metamist then has no CRAM analysis for the SG, so the
sample sheet asks for it to be aligned again. Under [ADR 0009](0009-never-overwrite-preflight-and-publishdir-no-overwrite.md)
the pre-flight check stopped the whole run on that CRAM until someone registered or deleted it by hand.

Now, when both the CRAM and `.crai` are at the SG's destination, the pre-flight check takes the SG as
**unregistered**: its rows become one row of kind `unregistered` that points at that CRAM. It isn't
realigned. Only its `.somalier` is checked, and the run fails if that exists. The CRAM is then fingerprinted
and registered as one analysis, with the fingerprint as a secondary file
([ADR 0012](0012-one-analysis-per-cram.md)).
The CRAM is assumed to have been written against `ref_fasta`, because only nf-align publishes to that path.

This relaxes ADR 0009's rule that a run never registers a file it didn't write. Any other partial set of
outputs, such as a `.crai` without its CRAM, still fails the run.

## Considered Options

- **Keep failing the run**: one stray CRAM held up every other SG in the run.
- **Fingerprint, but don't register the CRAM**: the next run still sees the SG as unaligned, and fails on the
  `.somalier` this run published.

## Consequences

- If `ref_fasta` has changed since the CRAM was written, somalier can't decode the CRAM, and the SG is
  reported in `failed.csv` at `somalier`. The CRAM stays unregistered.
- The CRAM isn't checked any further. A truncated CRAM from a failed copy would be registered, so long as its
  `.crai` was there too.
- `crams.csv` lists only the CRAMs this run aligned, not the ones it found.
- `failed.csv` reports the SG once, at `somalier` or `register`.
