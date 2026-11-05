# A run only registers files it wrote

nf-align registers a file in Metamist only when this run (or an earlier attempt of the same session, on
`-resume`) published it. An output already at an SG's destination is a clash, and the pre-flight check fails
the run on it. That includes a CRAM and `.crai` an earlier run published but never registered.

This supersedes [ADR 0011](0011-register-and-fingerprint-unregistered-crams.md), which took such a CRAM as
**unregistered**, fingerprinted it, and registered it without realigning. That had the run vouch for a file
it had never seen written: it couldn't tell a complete CRAM from a truncated copy, and had to assume which
reference the CRAM was written against. It restores
[ADR 0009](0009-never-overwrite-preflight-and-publishdir-no-overwrite.md)'s original rule.

## Considered Options

- **Register unregistered CRAMs (ADR 0011)**: no manual step, but it registers files whose provenance the
  run can't check.

## Consequences

- An SG that failed at `somalier` or `register` leaves a published, unregistered CRAM. The next run for it
  stops at the pre-flight check, until someone registers or deletes the files by hand.
- One stray file holds up every SG in the run, as it did before ADR 0011.
- ADR 0009's accepted gap remains: two runs launched close together for the same SG can both pass the
  check, and the second then registers the first run's CRAM.
