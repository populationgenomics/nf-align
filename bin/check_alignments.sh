#!/usr/bin/env bash
# Fail unless FILE is a complete BAM/CRAM with at least one record. An OOM-killed aligner can
# leave a well-formed but header-only file behind, which quickcheck alone would accept.
# usage: check_alignments.sh FILE [REF_FASTA]
set -euo pipefail

file=$1
ref=${2:-}

samtools quickcheck "$file"
# only the first record is decoded; `|| true` absorbs the SIGPIPE when head exits
first=$(samtools view ${ref:+-T "$ref"} "$file" | head -n1 || true)
if [ -z "$first" ]; then
  echo "check_alignments.sh: $file has no alignment records" >&2
  exit 1
fi
